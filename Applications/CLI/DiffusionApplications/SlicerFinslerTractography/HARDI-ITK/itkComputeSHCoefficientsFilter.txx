/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeSHCoefficientsFilter.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkComputeSHCoefficientsFilter_txx
#define _itkComputeSHCoefficientsFilter_txx

#include "itkComputeSHCoefficientsFilter.h"
#include "itkImageRegionIterator.h"
#include "math.h"
#include "vnl/algo/vnl_matrix_inverse.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
ComputeSHCoefficientsFilter< TInputImage, TOutputImage >
::ComputeSHCoefficientsFilter()
{
   if( TInputImage::ImageDimension != 3 )
      itkExceptionMacro( << "This filter is only supported for input dimension 3" );
   m_List.clear();
   m_BValue = 1000.0f;
   m_Lambda = 0.006f;
   m_DecompositionType = cOPDT;
   m_L = 6;;
   m_Baselines.clear();
   m_Gradients.clear();
}
   
template< class TInputImage, class TOutputImage >
void ComputeSHCoefficientsFilter< TInputImage, TOutputImage >
::GenerateOutputInformation( void )
{
   Superclass::GenerateOutputInformation();
   this->GetOutput()->SetVectorLength( shmaths::getNumberOfEvenAssociatedLegendrePolynomials( this->GetL() ) );
}
   
template< class TInputImage, class TOutputImage >
void ComputeSHCoefficientsFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData( void )
{
   Superclass::BeforeThreadedGenerateData();
   // We have to configure and precompute all the matrixes involved in the computation
   // of SH coefficients. Note that this is done only once, so the performance of the
   // algorithm is not compromised.
   
   // First, determine which channels are actually gradients:
   unsigned int total = m_List.size();
   m_Baselines.clear();
   m_Gradients.clear();
   for( unsigned int k=0; k<total; ++k ){
      double norm = ::sqrt( m_List[k][0]*m_List[k][0] + m_List[k][1]*m_List[k][1] + m_List[k][2]*m_List[k][2] );
      if( norm>1e-4 )
         m_Gradients.push_back( k );
      else
         m_Baselines.push_back( k );
   }
   // So the actual number of gradient directions is:
   unsigned int N = m_Gradients.size();
   // Now, generate gx, gy, and gz for each gradient direction:
   double* gx = new double[N];
   double* gy = new double[N];
   double* gz = new double[N];
   for( unsigned int k=0; k<N; ++k ){
      gx[k] = m_List[ m_Gradients[k] ][0];
      gy[k] = m_List[ m_Gradients[k] ][1];
      gz[k] = m_List[ m_Gradients[k] ][2];
   }
   // Now, compute the spherical coordinates corresponding to these gradient directions:
   double* theta = new double[N];
   double* phi   = new double[N];
   shmaths::computeShericalCoordsFromCartesian( gx, gy, gz, theta, phi, N );
   // Cartesian coordinates are no longer necessary:
   delete[] gx;
   delete[] gy;
   delete[] gz;
   // From the spherical coordinates, we can compute all required matrixes (standard convention):
   MatrixType B;
   MatrixType L;
   shmaths::computeSHMatrixSymmetric( N, theta, phi, this->GetL(), B  ); // B has size N x vectorLengthOfTheOutput
   shmaths::computeSHEigMatrixSymmetric( this->GetL(), L );              // L has size vectorLengthOfTheOutput x vectorLengthOfTheOutput
   // Generate the matrix solving the regularized LS problem:
   m_LS  = ( B.transpose() * B + (L*L)*m_Lambda );
   // Note that the inversion of m_LS is done only once, so performance is secundary.
   // In particular, we can use the SVD-based vnl_matrix_inverse fnction instead of
   // vnl_inverse
   m_LS  = vnl_matrix_inverse<double>( m_LS );
   // And, finally:
   m_LS *= ( B.transpose() );
   // If only the ADC is needed, it is enough to compute this matrix; however, we
   // consider the case when either Q-Balls, the cOPDT or the pOPDT need to be
   // computed. In these cases, it is ncessary to compute the FRT matrix for SH:
   if( m_DecompositionType != ADC   &&   m_DecompositionType != ATT ){
      MatrixType F;
      shmaths::computeSHFRTMatrixSymmetric( this->GetL(), F );
      if( m_DecompositionType == QBall )
         m_LS = F*m_LS;
      else
         m_LS = ( F*(L*m_LS) )*( 0.0625/PI/PI );
   }
   // Delete the memory previously allocated:
   delete[] theta;
   delete[] phi;
}
   
template< class TInputImage, class TOutputImage >
void ComputeSHCoefficientsFilter< TInputImage, TOutputImage >
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
   unsigned int R  = shmaths::getNumberOfEvenAssociatedLegendrePolynomials( this->GetL() );
   unsigned int N  = m_Gradients.size();
   unsigned int B  = m_Baselines.size();
   double inverseB = 1.0f/m_BValue;
   // Auxiliar matrices to use with the vnl matrix multiplication methods:
   MatrixType vi( N, 1 );
   MatrixType vo;
   
   for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it ){
      InputPixelType  ip = bit.Get();
      OutputPixelType op;
      op.SetSize( R );
      // Average all available baselines:
      double baseline = 0.0f;
      for( unsigned int k=0; k<B; ++k )
         baseline += ip[ m_Baselines[k] ];
      // Check for nans:
      if( baseline>1e-6 )
         baseline        = B/baseline;
      else{
         // This is a background point, so attenuation hould be close to 0
         // Hence, we put a very high value of the baseline to force such
         // attenuation. Note that "baseline" is in fact the inverse of
         // of S_0, so we use: 1/S_0 = 1/1e6 = 1e-6
         baseline        = 1e-6;
      }
      // Compute the auxiliar vector to multiply:
      for( unsigned int k=0; k<N; ++k ){
         // Normalize to compute the attenuation signal, E(q) = S_i/S_0:
         double lval = ip[ m_Gradients[k] ] * baseline;
         // Depending on the type of SH expansion to compute, we need
         // to perform different kinds of operations on E(q).
         if( lval<1e-6 )  // Avoid numerical issues
            lval = 1e-6;
         if( lval>1 )     // Avoid numerical issues
            lval = 1;
         switch( m_DecompositionType ){
            case ADC:
               vi( k, 0 ) = -::log( lval ) * inverseB;
               break;
            case QBall:
               vi( k, 0 ) = lval;
               break;
            case cOPDT:
               vi( k, 0 ) = -shmaths::Ein(  -::log( lval )  );
               break;
            case pOPDT:
               lval = -::log( lval );
               if ( lval<1e-6 )
                  lval = 1e-6;
               vi( k, 0 ) = ::log( lval );
               break;
            case ATT:
               vi( k, 0 ) = lval;
               break;
            default:
               break;
         }
      }
      // Multiply by the LS matrix to perform LS fitting:
      vo = m_LS * vi;
      // If we are computing either the cOPDT or the pOPDT, it is necessary to correct
      // the first element of vi, corresponding to the integral of the radial part of
      // the Laplacian (or, alternatively, the normalization constant so that the 
      // orientation function has integral 1 and hence it is a true PDF).
      if( m_DecompositionType==cOPDT || m_DecompositionType==pOPDT )
         vo( 0, 0 ) = 0.25f/PI;
      // Place these values in the output pixel:
      for( unsigned int k=0; k<R; ++k )
         op[k] = (float)( vo[k][0] );
      it.Set( op );
   }   
}
   

   
} // end namespace itk


#endif
