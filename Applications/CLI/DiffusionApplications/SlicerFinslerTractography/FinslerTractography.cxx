/*=========================================================================

  Program:   Diffusion Applications
  Language:  C++
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/dwiNoiseFilter/dwiNoiseFilter.cxx $
  Date:      $Date: 2008-11-25 18:46:58 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7972 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#include <iostream>
#include <algorithm>
#include <string>
#include <itkMetaDataObject.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkVectorImage.h>

#ifdef _WIN32
// to pick up M_SQRT2 and other nice things...
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "itkPluginUtilities.h"
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include "itkCastImageFilter.h"
#include "FinslerTractographyCLP.h"
//--------------------------------------------------
// Specific includes:
#include "itkParallelfastSweeping.h"
#include "HARDI-ITK/sphericalHarmonics.h"
#include "HARDI-ITK/itkComputeSHCoefficientsFilter.h"
#include "HARDI-ITK/itkComputeLocalCostFromSH.h"
//--------------------------------------------------
#define DIMENSION 3

template<class PixelType>
int DoIt( int argc, const char * argv[], PixelType )
{
   PARSE_ARGS;
   // do the typedefs
   typedef itk::VectorImage<PixelType,DIMENSION>                       DiffusionImageType;
   typedef itk::Image<float,DIMENSION>                                 ScalarImageType;
   typedef itk::Image<unsigned char,DIMENSION>                         LabelImageType;
   typedef itk::Image<itk::CovariantVector<float,DIMENSION>,DIMENSION> VectorImageType;
   
   typedef itk::CovariantVector<double,DIMENSION>       CovariantVectorType;
   
   typedef itk::ImageFileReader<LabelImageType>     LabelReaderType;
   
   typedef itk::VectorImage<float,DIMENSION>                                SHImageType;      // Each voxel is the list of SH coefficients
   typedef itk::VectorImage<float,DIMENSION>                                CostsImageType;   // Each voxel is the directional cost
   // Filter to compute the SH coefficients from the DWI:
   typedef itk::ComputeSHCoefficientsFilter<DiffusionImageType,SHImageType> SHComputeType;
   typedef typename SHComputeType::Pointer                                  SHComputePointer;
   typedef typename SHComputeType::GradientType                             GradientType;
   typedef typename SHComputeType::ListOfGradientsType                      ListOfGradientsType;
   // Filter to compute the local costs from the SH coefficients
   typedef itk::ComputeLocalCostFromSH<SHImageType,CostsImageType>          CostsComputeType;
   typedef typename CostsComputeType::Pointer                               CostsComputePointer;
   
   typedef itk::ParallelFastSweeping<CostsImageType,ScalarImageType,LabelImageType> FinslerDistanceComputeType;
   typedef typename FinslerDistanceComputeType::Pointer                             FinslerDistanceComputePointer;
   
   typedef itk::ImageFileWriter<ScalarImageType> ScalarWriterType;
   typedef itk::ImageFileWriter<VectorImageType> VectorWriterType;
   
   //================================================================================
   //================================================================================
   /** I- READ THE INPUT DWI VOLUME AND INTERPRET THE METADATA */
   
   std::vector< CovariantVectorType >                   diffusionDirections;
   
   typedef itk::ImageFileReader<DiffusionImageType>     FileReaderType;
   typename FileReaderType::Pointer reader = FileReaderType::New();
   reader->SetFileName( inputDWI.c_str() );
   try{
      reader->Update();
   }
   catch( itk::ExceptionObject& e ){
      std::cerr << "Cannot read seeds volume" << std::endl;
      return EXIT_FAILURE;
   }
   
   itk::MetaDataDictionary            imgMetaDictionary = reader->GetMetaDataDictionary();    
   std::vector<std::string>           imgMetaKeys       = imgMetaDictionary.GetKeys();
   std::vector<std::string>::iterator itKey             = imgMetaKeys.begin();
   std::string                        metaString;
   
   typedef itk::MetaDataDictionary DictionaryType;
   const DictionaryType & dictionary = reader->GetMetaDataDictionary();
   
   typedef itk::MetaDataObject< std::string > MetaDataStringType;
   DictionaryType::ConstIterator itr = dictionary.Begin();
   DictionaryType::ConstIterator end = dictionary.End();
   
   double       dBValue      = 1000;
   int          iFoundBValue = 0;
   unsigned int channels     = 0;
   while( itr != end ){
      itk::MetaDataObjectBase::Pointer entry = itr->second;
      MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;
      if( entryvalue ){
         ::size_t pos = itr->first.find("DWMRI_gradient");
         if ( pos != std::string::npos ){
            std::string tagkey = itr->first;
            std::string tagvalue = entryvalue->GetMetaDataObjectValue();
            double dx[DIMENSION];
            std::sscanf(tagvalue.c_str(), "%lf %lf %lf\n", &dx[0], &dx[1], &dx[2]);
            CovariantVectorType dGrad = (CovariantVectorType)(dx);
            if( dGrad.GetNorm()>1e-6 )
               dGrad.Normalize();
            diffusionDirections.push_back( dGrad );
            ++channels;
         }
         else{
            // try to find the b-value
            ::size_t posB = itr->first.find("DWMRI_b-value");
            if ( posB != std::string::npos ){
               std::string tagvalue = entryvalue->GetMetaDataObjectValue();
               std::sscanf(tagvalue.c_str(), "%lf\n", &dBValue );
               iFoundBValue = 1;
            }
         }
      }
      ++itr;
   }
   // find the first zero baseline image and use it for the noise estimation
   ::size_t iNrOfDWIs = diffusionDirections.size();
   ::size_t iFirstBaseline = std::string::npos;
   for ( ::size_t iI=0; iI<iNrOfDWIs; iI++ ){
      if ( diffusionDirections[iI].GetNorm()==0 ){
         iFirstBaseline = iI;
         break;
      }
   }
   if ( iFirstBaseline == std::string::npos ){
      std::cout << "Did not find an explicit baseline image." << std::endl;
      std::cout << "Treating the first volume as the baseline volume." << std::endl;
      iFirstBaseline = 0;
   }
   //================================================================================
   //================================================================================
   /** II- READ THE VOLUME OF SEEDING POINTS */
   typename LabelReaderType::Pointer seedsReader = LabelReaderType::New();
   seedsReader->SetFileName( inputSeeds.c_str() );
   try{
      seedsReader->Update();
      // Check if this has the appropriate size
      if( seedsReader->GetOutput()->GetLargestPossibleRegion() != reader->GetOutput()->GetLargestPossibleRegion() ){
         std::cout << "The seeds volume has not an adequare size. I cannot proceed" << std::endl;
         return EXIT_FAILURE;
      }
   }
   catch( itk::ExceptionObject& e ){
      std::cerr << "Cannot read seeds volume" << std::endl;
      return EXIT_FAILURE;
   }
   //================================================================================
   //================================================================================
   /** III- READ THE VOLUME WITH THE MASK */
   typename LabelReaderType::Pointer maskReader = LabelReaderType::New();
   maskReader->SetFileName( inputMask.c_str() );
   bool useMask;
   if( inputMask.length() > 0 ){
      try{
         maskReader->Update();
         // Check if this has the appropriate size
         if( maskReader->GetOutput()->GetLargestPossibleRegion() != reader->GetOutput()->GetLargestPossibleRegion() ){
            std::cout << "The mask has a weird size. I'm not using it..." << std::endl;
            useMask = false;
         }
         else
            useMask = true;
      }
      catch( itk::ExceptionObject& e ){
         std::cout << "I dindn't find a mask. Running without it..." << std::endl;
         useMask = false;
      }
   }
   else{
      std::cout << "I'm not using the mask" << std::endl;
      useMask = false;
   }
   //================================================================================
   //================================================================================
   /** IV- COMPUTE THE VOLUME OF SH AND THE DIRECTIONAL COSTS */
   // Create filters and do the job:
   SHComputePointer              shCompute      = SHComputeType::New();
   CostsComputePointer           costsCompute   = CostsComputeType::New();
   // SET THE PARAMETERS AND THE INPUT TO THE FILTER COMPUTING THE SH COEFFICIENTS
   // Set the list of gradient directions (retrieve from the DWI volume):
   shCompute->SetList( diffusionDirections );
   // Set the b-value of the acquisition (retireve from the DWI volume):
   shCompute->SetBValue( dBValue );
   // Set the "lambda" regularization parameter (retrieve from the GUI):
   shCompute->SetLambda( iLambda );
   // Set the degree of the SH expansions (retireve from the GUI;
   //       future work: check or automatically determine from the number of DWI channels):
   shCompute->SetL( iLSH );
   // Set the DWI-related measure (right now, always the attenaution signal
   //       future work: allow different ways to compute the local cost):
   shCompute->SetComputeATT();
   // Set the input:
   shCompute->SetInput( reader->GetOutput() );
   //-----------------------------------------------------------------------------------------------------------------------------------
   // SET THE PARAMETERS AND THE INPUT TO THE FILTER COMPUTING THE LOCAL COST
   // Default neighboring directions (Future work: accept these directions as a parameter):
   typedef FinslerDistanceComputeType::DirectionType NeighborDirectionType;
   std::vector<NeighborDirectionType> ndirs;
   NeighborDirectionType              ndir;
   for( int z=-1; z<=1; ++z ){
      for( int y=-1; y<=1; ++y ){
         for( int x=-1; x<=1; ++x ){
            if( x!=0 || y!=0 || z!=0 ){
               ndir[0] = x;
               ndir[1] = y;
               ndir[2] = z;
               ndir.Normalize();
               ndirs.push_back( ndir );
            }
         }
      }
   }
   // Set the neighboring directions:
   costsCompute->SetNeighboringDirections( ndirs );
   // Set the input:
   costsCompute->SetInput( shCompute->GetOutput() );
   //================================================================================
   //================================================================================
   /** V- DO FAST SWEEPING */
   FinslerDistanceComputePointer fastSweeping = FinslerDistanceComputeType::New();
   fastSweeping->SetInput( costsCompute->GetOutput() );
   //-----------------------------------------------------------------------------------------------------------------------------------
   // SET THE PARAMETERS AND THE INPUT TO THE FILTER PERFORMING FAST SWEEPING
   // Set the maximum number of iterations (retrieve from the GUI):
   fastSweeping->SetMaxIters( iMaxIters );
   // Set the cost factor for the stop criterion (retrieve log-value from the GUI):
   double factor = (double)iCostFraction;
   factor = ::exp( factor * ::log(10.0f) );
   fastSweeping->SetCostFactor( factor );
   // Set the neighboring directions:
   fastSweeping->SetNeighboringDirections( ndirs );
   // Set the input
   fastSweeping->SetInput( costsCompute->GetOutput() );
   // Set the seeding points:
   std::vector<ScalarImageType::IndexType > seeds;
   seeds.clear();
   typedef itk::ImageRegionConstIteratorWithIndex<LabelImageType> IteratorType;
   IteratorType it( seedsReader->GetOutput(), seedsReader->GetOutput()->GetLargestPossibleRegion() );
   for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
      if( it.Get() )
         seeds.push_back( it.GetIndex() );
   }
   // Make sure at least one seed is present:
   if( seeds.size()==0 ){
      std::cerr << "I have no seeds to work with!" << std::endl;
      return EXIT_FAILURE;
   }
   // Set the seeding points of the Finsler filter:
   fastSweeping->SetSeedingPoints( seeds );
   //-----------------------------------------------------------------------------------------------------------------------------------
   // Set the mask
   if( useMask )
      fastSweeping->SetMask( maskReader->GetOutput() );
   //-----------------------------------------------------------------------------------------------------------------------------------
   // Prepare to write
   typename ScalarWriterType::Pointer writer = ScalarWriterType::New();
   writer->SetFileName( outputCost.c_str() );
   writer->SetInput( fastSweeping->GetOutput() );
   try{
      writer->Update();
   }
   catch( itk::ExceptionObject& e ){
      std::cerr << "Cannot update the filter, something went wrong..." << std::endl;
      return EXIT_FAILURE;
   }
   //================================================================================
   //================================================================================
   /** VI- WRITE OUT THE MAP OF ARRIVAL DIRECTIONS */
   if( outputDirections.length() > 0 ){
      typename VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
      vectorWriter->SetFileName( outputDirections.c_str() );
      vectorWriter->SetInput( fastSweeping->GetOptimumDirectionsMap() );
      try{
         vectorWriter->Update();
      }
      catch( itk::ExceptionObject& e ){
         std::cout << "I cannot write the map of optimal output directions" << std::endl;
      }
   }
   else{
      std::cout << "I am not writting the map of optimum directions" << std::endl;
   }

   //================================================================================
   //================================================================================
   return EXIT_SUCCESS;
}

int main( int argc, const char * argv[] )
{
  
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  //try
  // {
  itk::GetImageType (inputDWI, pixelType, componentType);

  // This filter handles all types
    
  switch (componentType)
    {
#ifndef WIN32
    case itk::ImageIOBase::UCHAR:
      return DoIt( argc, argv, static_cast<unsigned char>(0));
      break;
    case itk::ImageIOBase::CHAR:
      return DoIt( argc, argv, static_cast<char>(0));
      break;
#endif
    case itk::ImageIOBase::USHORT:
      return DoIt( argc, argv, static_cast<unsigned short>(0));
      break;
    case itk::ImageIOBase::SHORT:
      return DoIt( argc, argv, static_cast<short>(0));
      break;
    case itk::ImageIOBase::UINT:
      return DoIt( argc, argv, static_cast<unsigned int>(0));
      break;
    case itk::ImageIOBase::INT:
      return DoIt( argc, argv, static_cast<int>(0));
      break;
#ifndef WIN32
    case itk::ImageIOBase::ULONG:
      return DoIt( argc, argv, static_cast<unsigned long>(0));
      break;
    case itk::ImageIOBase::LONG:
      return DoIt( argc, argv, static_cast<long>(0));
      break;
#endif
    case itk::ImageIOBase::FLOAT:
      std::cout << "FLOAT type not currently supported." << std::endl;
      break;
    case itk::ImageIOBase::DOUBLE:
      std::cout << "DOUBLE type not currently supported." << std::endl;
      break;
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    default:
      std::cout << "unknown component type" << std::endl;
      break;
    }

  return EXIT_SUCCESS;
}


