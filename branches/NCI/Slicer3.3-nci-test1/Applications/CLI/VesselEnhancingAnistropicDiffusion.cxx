/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/GradientAnisotropicDiffusion.cxx $
  Language:  C++
  Date:      $Date: 2006-12-18 09:29:35 -0500 (Mon, 18 Dec 2006) $
  Version:   $Revision: 1855 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"

#include "itkAnisotropicDiffusionVesselEnhancementImageFilter.h"

#include "itkPluginUtilities.h"
#include "VesselEnhancingAnistropicDiffusionCLP.h"

template<class T> int DoIt( int argc, char * argv[], T )
{

  PARSE_ARGS;

  //Pixel Types
  typedef    double   InputPixelType;
  typedef    double   OutputPixelType;
//  typedef    T       OutputPixelType;

  //Image Types
  typedef itk::Image< InputPixelType,  3 >   InputImageType;
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

  //Read Image
  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  itk::PluginFilterWatcher watchReader(reader, "Read Volume",
                                       CLPProcessInformation);
  reader->SetFileName( inputVolume.c_str() );

  //Filter Image
  typedef itk::AnisotropicDiffusionVesselEnhancementImageFilter<
               InputImageType, InputImageType >  FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  itk::PluginFilterWatcher watchFilter(filter, "Vessel Enhancing Anisotropic Diffusion",
    CLPProcessInformation);

  filter->SetInput( reader->GetOutput() );   
  filter->SetAlpha( Alpha );
  filter->SetBeta( Beta );
  filter->SetGamma( Gamma );
  filter->SetC( C );
  filter->SetScaleVesselnessMeasure( ScaleVesselnessMeasure );
  filter->SetSigmaMin( StartScale );
  filter->SetSigmaMax( StopScale );
  filter->SetNumberOfSigmaSteps( NumberOfScales ); 
  filter->SetEpsilon( Epsilon); 
  filter->SetWStrength( Omega ); 
  filter->SetSensitivity( Sensitivity ); 
  filter->SetTimeStep( TimeStep );
  filter->SetNumberOfIterations( NumberOfIterations );

  //Write Image
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  itk::PluginFilterWatcher watchWriter(writer, "Write Volume",
                                   CLPProcessInformation);
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{

  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType (inputVolume, pixelType, componentType);

    // This filter handles all types

    switch (componentType)
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0));
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0));
        break;
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
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0));
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0));
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0));
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0));
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }
  catch( itk::ExceptionObject &excep)
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

