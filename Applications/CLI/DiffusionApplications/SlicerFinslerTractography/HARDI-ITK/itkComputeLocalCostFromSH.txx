/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeLocalCostFromSH.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkComputeLocalCostFromSH_txx
#define _itkComputeLocalCostFromSH_txx

#include "itkComputeLocalCostFromSH.h"
#include "itkImageRegionIterator.h"
#include "math.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
ComputeLocalCostFromSH< TInputImage, TOutputImage >
::ComputeLocalCostFromSH()
{
   if( TInputImage::ImageDimension != 3 )
      itkExceptionMacro( << "This filter is only supported for input dimension 3" );
   m_NeighboringDirections.resize(0);
}
   
template< class TInputImage, class TOutputImage >
void ComputeLocalCostFromSH< TInputImage, TOutputImage >
::GenerateOutputInformation( void )
{
   Superclass::GenerateOutputInformation();
   this->GetOutput()->SetVectorLength( m_NeighboringDirections.size() );
}
      
template< class TInputImage, class TOutputImage >
void ComputeLocalCostFromSH< TInputImage, TOutputImage >
::BeforeThreadedGenerateData( void )
{
   // Call superclass' implementation for this method:
   Superclass::BeforeThreadedGenerateData();
   // Compute a couple of parameters of interest and check:
   unsigned int N = m_NeighboringDirections.size();
   if( N==0 )
      itkExceptionMacro( << "You must set the directions to compute the costs before executing the filter" );
   unsigned int M = this->GetInput()->GetVectorLength();
   // From M, we need to compute the order of the SH:
   unsigned int L = 0;
   switch(M){
      case 1:
         L = 0;
         break;
      case 6:
         L = 2;
         break;
      case 15:
         L = 4;
         break;
      case 28:
         L = 6;
         break;
      case 45:
         L = 8;
         break;
      default:
         itkExceptionMacro( << "Only degrees 0, 2, 4, 6, and 8 are supported. (There are " << M << " coefficients)" );
         break;
   }
   
   // The directions to compute the costs are given in cartesian coordinates. We
   // need to use spherical coordinates in order to compute the matrix to convert
   // from SH coefficients to actual costs. This is done by the shmaths code:
   double* gx = new double[N]; // Old-style c vector
   double* gy = new double[N]; // Old-style c vector
   double* gz = new double[N]; // Old-style c vector
   // Retrieve the cartesian coordinates in standard c vectors:
   for( unsigned int k=0; k<N; ++k ){
      if( m_NeighboringDirections[k].GetNorm()>1e-6 )
         m_NeighboringDirections[k].Normalize();
      gx[k] = m_NeighboringDirections[k][0];
      gy[k] = m_NeighboringDirections[k][1];
      gz[k] = m_NeighboringDirections[k][2];
   }
   // Now, compute the spherical coordinates corresponding to these gradient directions:
   double* theta = new double[N];
   double* phi   = new double[N];
   shmaths::computeShericalCoordsFromCartesian( gx, gy, gz, theta, phi, N );
   // Cartesian coordinates are no longer necessary:
   delete[] gx; delete[] gy; delete[] gz;
   // From the spherical coordinates, we can compute all required matrixes (standard convention):
   shmaths::computeSHMatrixSymmetric( N, theta, phi, L, m_B  ); // B has size N x M
   shmaths::computeSHFRTMatrixSymmetric( L, m_FRT );            // FRT matrix has size M x M
   // The matrix m_B is the one required to compute the costs for each direction from
   // the SH coefficients at each voxel. Delete the memory previously allocated:
   delete[] theta;
   delete[] phi;
}
   
template< class TInputImage, class TOutputImage >
void ComputeLocalCostFromSH< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
   // Iterators:
   ImageRegionConstIterator<InputImageType>        bit;  // Iterator for the input image
   ImageRegionIterator<OutputImageType>            it;   // Iterator for the output image
   // Input and output
   InputImageConstPointer   input   =  this->GetInput();
   OutputImagePointer       output  =  this->GetOutput();
   
   bit = ImageRegionConstIterator<InputImageType>( input,  outputRegionForThread );
   it  = ImageRegionIterator<OutputImageType>(     output, outputRegionForThread );
   
   // Precompute a couple of parameters of interest:
   unsigned int N = m_NeighboringDirections.size();
   unsigned int M = this->GetInput()->GetVectorLength();
   // Auxiliar value to perform matrix products:
   MatrixType pixel(M,1);
   OutputPixelType op(N);
   for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it ){
      // Get the pixel and convert to a suitable format:
      InputPixelType  ip = bit.Get();
      for( unsigned int m=0; m<M; ++m )
         pixel[m][0] = ip[m];
      
      //OLD IMPLEMENTATION, CLONING JOHN'S CODE
      // Convert to costs:
      MatrixType att = m_B*pixel;
      MatrixType nor = m_B*m_FRT*pixel;
      // Compute the final output and set:
      for( unsigned int n=0; n<N; ++n ){
         double num   = ( att[n][0]>0 ? att[n][0] :        0       );
         double den   = ( nor[n][0]>0 ? nor[n][0] : 0.01*att[n][0] );
         double value = ( num + 0.00001f ) / ( den + 0.00001f );
         op[n]        = (value*value*value);
      }
       
      /*
      // DUMB TEST FOR DEBUG
      double th = 0.5;
      for( unsigned int n=0; n<N; ++n ){
         DirectionType d1 = m_NeighboringDirections[n];
         d1.Normalize();
         DirectionType d2;
         d2[0] = 0.0f;
         d2[1] = 1.0f;
         d2[2] = 0.0f;
         double temp = d1*d2;
         temp  = 1.0f - temp;
         temp = ( temp>th ? temp : th );
         op[n] = ( temp * temp * temp );
      }
      */
      
      it.Set( op );
   }
}
   
} // end namespace itk


#endif
