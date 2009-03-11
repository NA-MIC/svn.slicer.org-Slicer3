/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleHessianLine3DFilter.txx,v $
  Language:  C++
  Date:      $Date:   $
  Version:   $Revision:   $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiScaleHessianLine3DFilter_txx
#define __itkMultiScaleHessianLine3DFilter_txx



#include "itkMultiScaleLine3DFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

#define EPSILON  1e-03

namespace itk
{

/**
 * Constructor
 */
template <typename TInputImage, typename TOutputImage >
MultiScaleLine3DFilter<TInputImage,TOutputImage>
::MultiScaleLine3DFilter()
{
  m_SigmaMin = 0.2;
  m_SigmaMax = 2.0;

  m_NumberOfSigmaSteps = 10;

  m_HessianFilter                = HessianFilterType::New();
  m_VesselnessFilter             = VesselnessFilterType::New();

  //Instantiate Update buffer
  m_UpdateBuffer                 = UpdateBufferType::New();
}

template <typename TInputImage, typename TOutputImage >
void
MultiScaleLine3DFilter<TInputImage,TOutputImage>
::AllocateUpdateBuffer()
{
  /* The update buffer looks just like the output and holds the best response
     in the  vesselness measure */
  
  typename TOutputImage::Pointer output = this->GetOutput();

  m_UpdateBuffer->SetSpacing(output->GetSpacing());
  m_UpdateBuffer->SetOrigin(output->GetOrigin());
  m_UpdateBuffer->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  m_UpdateBuffer->SetRequestedRegion(output->GetRequestedRegion());
  m_UpdateBuffer->SetBufferedRegion(output->GetBufferedRegion());
  m_UpdateBuffer->Allocate();
}


template <typename TInputImage, typename TOutputImage >
void
MultiScaleLine3DFilter<TInputImage,TOutputImage>
::GenerateData()
{
  // Allocate the output
  this->GetOutput()->SetBufferedRegion( 
                 this->GetOutput()->GetRequestedRegion() );
  this->GetOutput()->Allocate();

  // Allocate the buffer
  AllocateUpdateBuffer();

  typename InputImageType::ConstPointer input = this->GetInput();
 
  this->m_HessianFilter->SetInput( input );

  this->m_HessianFilter->SetNormalizeAcrossScale( false );
 
  double sigma = m_SigmaMin;

  int scaleLevel = 1;

  while ( sigma <= m_SigmaMax )
    {
    std::cout << "Computing vesselness for scale with sigma= " 
              << sigma << std::endl;

    m_HessianFilter->SetSigma( sigma );

    m_VesselnessFilter->SetInput ( m_HessianFilter->GetOutput() ); 
    m_VesselnessFilter->SetSigma( sigma );
    m_VesselnessFilter->Update();
 
    this->UpdateMaximumResponse();

    sigma  = this->ComputeSigmaValue( scaleLevel );

    scaleLevel++;
    } 

  //Write out the best response to the output image
  ImageRegionIterator<UpdateBufferType> 
               it(m_UpdateBuffer,m_UpdateBuffer->GetLargestPossibleRegion());
  it.GoToBegin();

  ImageRegionIterator<TOutputImage> oit(this->GetOutput(),
                          this->GetOutput()->GetLargestPossibleRegion());
  oit.GoToBegin();

  while(!oit.IsAtEnd())
    {
    oit.Value() = static_cast< OutputPixelType >( it.Get() );
    ++oit;
    ++it;
    }
}

template <typename TInputImage, typename TOutputImage >
void
MultiScaleLine3DFilter<TInputImage,TOutputImage>
::UpdateMaximumResponse()
{

  ImageRegionIterator<UpdateBufferType> 
            oit(m_UpdateBuffer,m_UpdateBuffer->GetLargestPossibleRegion());

  oit.GoToBegin();

  typedef typename VesselnessFilterType::OutputImageType   
                                         VesselnessOutputImageType;

  ImageRegionIterator<VesselnessOutputImageType> 
            it(m_VesselnessFilter->GetOutput(),
            this->m_VesselnessFilter->GetOutput()->GetLargestPossibleRegion());

  it.GoToBegin();

  while(!oit.IsAtEnd())
    {
    if( oit.Value() < it.Value() )
      {
      oit.Value() = it.Value();
      }
    ++oit;
    ++it;
    }
}

template <typename TInputImage, typename TOutputImage >
double
MultiScaleLine3DFilter<TInputImage,TOutputImage>
::ComputeSigmaValue( int ScaleLevel )
{
  double stepSize = 
     ( vcl_log( m_SigmaMax )  - vcl_log( m_SigmaMin) ) / m_NumberOfSigmaSteps;

  if( stepSize <= 1e-10 )
    {
    stepSize = 1e-10;
    } 

  return ( vcl_exp( vcl_log (m_SigmaMin) + stepSize * ScaleLevel) );
}

template <typename TInputImage, typename TOutputImage >
void
MultiScaleLine3DFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "SigmaMin:  " << m_SigmaMin << std::endl;
  os << indent << "SigmaMax:  " << m_SigmaMax  << std::endl;
}


} // end namespace itk
  
#endif
