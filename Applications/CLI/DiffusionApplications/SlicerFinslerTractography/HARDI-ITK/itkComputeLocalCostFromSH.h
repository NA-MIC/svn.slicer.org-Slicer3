/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeLocalCostFromSH.h,v $
  Language:  C++
  Date:      $Date: 2006/03/27 17:01:10 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComputeLocalCostFromSH_h
#define __itkComputeLocalCostFromSH_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "sphericalHarmonics.h"

namespace itk
{


template <class TInputImage, class TOutputImage>
class ITK_EXPORT ComputeLocalCostFromSH : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
   /** Standard class typedefs. */
   typedef ComputeLocalCostFromSH                Self;
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
   itkTypeMacro( ComputeLocalCostFromSH, ImageToImageFilter );
  
   /** Image typedef support. */
   typedef typename InputImageType::PixelType           InputPixelType;
   typedef typename OutputImageType::PixelType          OutputPixelType;
   typedef typename InputImageType::RegionType          InputImageRegionType;
   typedef typename InputImageType::SizeType            InputImageSizeType;
   typedef typename InputImageType::IndexType           InputImageIndexType;
   typedef typename OutputImageType::RegionType         OutputImageRegionType;
   
   /** Types to store the directions of interest */
   typedef itk::CovariantVector<float,TInputImage::ImageDimension> DirectionType;
   
   void SetNeighboringDirections( std::vector<DirectionType> dirs )
   {
      m_NeighboringDirections = dirs;
   }
   
   // Type for matrix computations (vnl matrix);
   typedef shmaths::SHMatrixType                        MatrixType;
protected:
   ComputeLocalCostFromSH();
   virtual ~ComputeLocalCostFromSH() {}
   void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
   void BeforeThreadedGenerateData( void );
   void GenerateOutputInformation( void );
private:
   ComputeLocalCostFromSH(const Self&);    // purposely not implemented
   void operator=(const Self&);            // purposely not implemented
   
   // The directions for which the cost will be computed:
   std::vector<DirectionType> m_NeighboringDirections;
   // The matrixes to compute the values from the SH coefficients:
   MatrixType                 m_B;   // SH matrix
   MatrixType                 m_FRT; // The Funk-Radon transform matrix
};


 
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeLocalCostFromSH.txx"
#endif

#endif
