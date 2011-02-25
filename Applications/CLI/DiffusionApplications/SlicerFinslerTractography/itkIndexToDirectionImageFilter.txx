/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkIndexToDirectionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkIndexToDirectionImageFilter_txx
#define _itkIndexToDirectionImageFilter_txx

#include "itkIndexToDirectionImageFilter.h"
#include "itkImageRegionIterator.h"
#include "math.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
IndexToDirectionImageFilter< TInputImage, TOutputImage >
::IndexToDirectionImageFilter()
{
   m_NeighboringDirections.resize(0);
}
      
template< class TInputImage, class TOutputImage >
void IndexToDirectionImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData( void )
{
   Superclass::BeforeThreadedGenerateData();
   if( m_NeighboringDirections.size()==0 )
      itkExceptionMacro( << "You must set the directions to use before updating" );
}
   
template< class TInputImage, class TOutputImage >
void IndexToDirectionImageFilter< TInputImage, TOutputImage >
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
   
   DirectionType zeroDisplacementValue;
   zeroDisplacementValue.Fill( 0.0f );
   // Precompute the number of directions:
   unsigned int N = m_NeighboringDirections.size();
   
   for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it ){
      unsigned int pos = static_cast<unsigned int>( bit.Get() );
      if( pos<N ){
         it.Set( static_cast<OutputPixelType>( m_NeighboringDirections[pos] ) );
      }
      else{
         it.Set( zeroDisplacementValue );
      }
   }   
}
   

   
} // end namespace itk


#endif
