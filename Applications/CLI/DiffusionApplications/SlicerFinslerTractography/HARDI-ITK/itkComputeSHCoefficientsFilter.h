/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeSHCoefficientsFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/03/27 17:01:10 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComputeSHCoefficientsFilter_h
#define __itkComputeSHCoefficientsFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include <vector>
#include "itkVector.h"
#include "sphericalHarmonics.h"
#include "itkCovariantVector.h"

namespace itk
{


template <class TInputImage, class TOutputImage>
class ITK_EXPORT ComputeSHCoefficientsFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
   /** Standard class typedefs. */
   typedef ComputeSHCoefficientsFilter                     Self;
   /** Convenient typedefs for simplifying declarations. */
   typedef TInputImage                           InputImageType;
   typedef typename InputImageType::Pointer      InputImagePointer;
   typedef typename InputImageType::ConstPointer InputImageConstPointer;
   typedef TOutputImage                          OutputImageType;
   typedef typename OutputImageType::Pointer     OutputImagePointer;
   
   /** Standard class typedefs. */
   typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
   typedef SmartPointer<Self>                                   Pointer;
   typedef SmartPointer<const Self>                             ConstPointer;

   /** Method for creation through the object factory. */
   itkNewMacro(Self);

   /** Run-time type information (and related methods). */
   itkTypeMacro( ComputeSHCoefficientsFilter, ImageToImageFilter );
  
   /** Image typedef support. */
   typedef typename InputImageType::PixelType           InputPixelType;
   typedef typename OutputImageType::PixelType          OutputPixelType;
   typedef typename InputImageType::RegionType          InputImageRegionType;
   typedef typename InputImageType::SizeType            InputImageSizeType;
   typedef typename InputImageType::IndexType           InputImageIndexType;
   typedef typename OutputImageType::RegionType         OutputImageRegionType;
   typedef typename OutputPixelType::ValueType          ScalarType;

   // Gradient type:
   typedef itk::CovariantVector<double,3>               GradientType;
   // Type to store all the gradient directions
   typedef std::vector< GradientType >                  ListOfGradientsType;
   // Type to store which indexes are baselines and which of them are gradients:
   typedef std::vector< unsigned int >                  IndicatorType;
   // Type for matrix computations (vnl matrix);
   typedef shmaths::SHMatrixType                        MatrixType;
   
   void SetList( ListOfGradientsType l ){
      m_List = l;
   }
   ListOfGradientsType GetList( void ){
      return m_List;
   }
   
   itkSetMacro( BValue, double );
   itkGetMacro( BValue, double );
   
   itkSetMacro( Lambda, double );
   itkGetMacro( Lambda, double );
   
   itkSetMacro( L, unsigned int );
   itkGetMacro( L, unsigned int );
   
   // This is to select which type of SH de composition we want to compute: ADC,
   // Q-Balls, cOPDT, or pOPDT:
   typedef enum{ ADC, QBall, cOPDT, pOPDT, ATT } DecompositionType;
   void SetComputeADC( void ){
      m_DecompositionType = ADC;
   }
   void SetComputeQBall( void ){
      m_DecompositionType = QBall;
   }
   void SetComputecOPDT( void ){
      m_DecompositionType = cOPDT;
   }
   void SetComputepOPDT( void ){
      m_DecompositionType = pOPDT;
   }
   void SetComputeATT( void ){
      m_DecompositionType = ATT;
   }
protected:
   ComputeSHCoefficientsFilter();
   virtual ~ComputeSHCoefficientsFilter() {}
   void GenerateOutputInformation( void );
   void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
   void BeforeThreadedGenerateData( void );
private:
   ComputeSHCoefficientsFilter(const Self&); // purposely not implemented
   void operator=(const Self&);    // purposely not implemented
   // The list of gradient directions:
   ListOfGradientsType     m_List;
   double                  m_BValue;
   double                  m_Lambda;
   unsigned int            m_L;
   IndicatorType           m_Baselines;
   IndicatorType           m_Gradients;
   MatrixType              m_LS;
   DecompositionType       m_DecompositionType;
};


 
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeSHCoefficientsFilter.txx"
#endif

#endif
