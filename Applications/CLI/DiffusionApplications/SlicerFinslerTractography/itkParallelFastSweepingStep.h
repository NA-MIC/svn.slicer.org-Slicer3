/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParallelFastSweepingStep.h,v $
  Language:  C++
  Date:      $Date: 2011-01-11 $
  Version:   $Revision: 1.0 $

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkParallelFastSweepingStep_h
#define __itkParallelFastSweepingStep_h


#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkArray.h"
#include "itkArray2D.h"
#include "vector.h"

namespace itk
{
/** \class ParallelFastSweepingStep
 * \brief Applies an iteration of the fast sweeping algorithm
 *
 */
template <class TInputImage, class TOutputImage, class TVectorImage, class TMaskImage = Image<unsigned char,TInputImage::ImageDimension> >
class ITK_EXPORT ParallelFastSweepingStep : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
   /** Extract dimension from input and output image. */
   itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
   itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);
   
   /** Convenient typedefs for simplifying declarations. */
   typedef TInputImage                       InputImageType;
   typedef typename InputImageType::Pointer  InputImagePointer;
   typedef TOutputImage                      OutputImageType;
   typedef typename OutputImageType::Pointer OutputImagePointer;
   
   /** Standard class typedefs. */
   typedef ParallelFastSweepingStep                             Self;
   typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
   typedef SmartPointer<Self>                                   Pointer;
   typedef SmartPointer<const Self>                             ConstPointer;
   
   /** Method for creation through the object factory. */
   itkNewMacro(Self);
   
   /** Run-time type information (and related methods). */
   itkTypeMacro(ParallelFastSweepingStep, ImageToImageFilter);
   
   /** Image typedef support. */
   typedef typename InputImageType::PixelType    InputPixelType;
   typedef typename OutputImageType::PixelType   OutputPixelType;
   
   typedef typename InputImageType::RegionType   InputImageRegionType;
   typedef typename OutputImageType::RegionType  OutputImageRegionType;
   
   typedef typename InputImageType::SizeType     InputSizeType;
   typedef typename InputImageType::IndexType    InputIndexType;
   
   typedef typename OutputImageType::SizeType    OutputSizeType;
   typedef typename OutputImageType::IndexType   OutputIndexType;
   typedef typename OutputImageType::OffsetType  OutputOffsetType;
   
   /** These types are used to describe the "arrival directions map" (see below) */
   typedef itk::Image<unsigned int,TInputImage::ImageDimension>          ArrivalDirectionsType;
   typedef typename ArrivalDirectionsType::Pointer                       ArrivalDirectionsPointer;
   /** These types are used to store the local costs for each arrival direction */
   typedef TVectorImage                                                  LocalCostType;
   typedef typename LocalCostType::ConstPointer                          LocalCostPointer;
   /** Type to store the mask for the input, if provided */
   typedef TMaskImage                                                    MaskType;
   typedef typename MaskType::Pointer                                    MaskPointer;
   
   void SetArrivalDirections( ArrivalDirectionsPointer arrivals ){
      m_ArrivalDirections = arrivals;
   }
   
   void SetLocalCost( LocalCostPointer costs ){
      m_LocalCost = costs;
   }
   
   void SetMask( MaskPointer mask ){
      m_Mask = mask;
   }
   
   void SetChosen( itk::Array<unsigned int>* chosen ){
      m_Chosen = chosen;
   }
   
   itkSetMacro( Neighbors, itk::Array2D<unsigned int> );
   itkGetMacro( Neighbors, itk::Array2D<unsigned int> );
   itkSetMacro( Weights,   itk::Array2D<float>        );
   itkGetMacro( Weights,   itk::Array2D<float>        );
   itkSetMacro( SplitDirection, unsigned int );
   itkGetMacro( SplitDirection, unsigned int );
   itkGetMacro( CostChange, OutputPixelType );

protected:
   ParallelFastSweepingStep();
   virtual ~ParallelFastSweepingStep() {}
   virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);
   void PrintSelf(std::ostream& os, Indent indent) const;
   void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId );
   void BeforeThreadedGenerateData();
   void AfterThreadedGenerateData();
   int SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion);
private:
   ParallelFastSweepingStep(const Self&); // purposely not implemented
   void operator=(const Self&);           // purposely not implemented
  
   itk::Array2D<unsigned int>      m_Neighbors;     // The neighbors to use for each discretized arriving direction
   itk::Array2D<float>             m_Weights;       // The interpolation weights to use for each arriving direction
   
   std::vector<OutputOffsetType>   m_Offsets;       // The set offsets for each neighbor
   OutputPixelType                 m_CostThreshold; // The threshold above which the cost is considered to be infinite
   
   /** This is a collection of vectors; each vector of the collection
    corresponds to a sweeping direction for the fast-sweeping algorithm;
    for each sweeping direction, the corresponding vector is a collection
    of indices indicating which sampled spatial directions are actually
    checked to update (or not) the cost at the current voxel. This way we aim
    to exploit the causaity of the solution to the differential equation to
    achieve a net speed-up. IMPORTANT: the user class is responsible to
    allocate and mantain the correspoding memory.*/
   itk::Array<unsigned int>*    m_Chosen;
   /** We keep a scalar image the same size as the input. Each pixel has an 
    unsigned integer value from 0 to the number of directions (i.e. the rows
    of m_Neighbors and m_Weights) minus 1. This value corresponds to the direction
    (or, more generally, the set of neighbors) producing the minimum accumulated
    cost at this voxel. For efficiency purposes, it is the responsibility of the
    user to allocate the required memory for this map and to fix the appropriate
    origin, spacing, regions size... Also, the user must initialize the values
    of the image, that will be iteratively updated by this filter.*/
   ArrivalDirectionsPointer        m_ArrivalDirections;
   /** This vector image has the same size as the input and the output. Its voxels
    are vectors with same number of components as discretized dirctions are
    considered, i.e., the same as the number of components of m_Neighbors. For
    a given component, the corresponding value represents the local cost at the
    current voxel associated to the direction being measured. Like the image right
    before, it is the responsibility of the user to allocate this image with the
    appropriate size and information, as well as initialize its values (they
    are not updated by this filter).*/
   LocalCostPointer                m_LocalCost;
   /** We can use a FA mask to accelerate computations and eliminate
    impossible or unlikely paths:*/
   MaskPointer                     m_Mask;
   /** This method overrides the default behavior for SplitRequestedRegion.
    Instead of splitting along the outermost available direction, the user
    may choose the preferred direction to split, which can be used to accelerate
    the propagation of the information from the seeding points and hence improve
    convergence.*/
   unsigned int                    m_SplitDirection;
   /** The change in the cost for the current iteration; we use a per-thread
    value and a global value to keep the filter thread-safe.*/
   itk::Array<OutputPixelType>     m_PerThreadCostChange;
   OutputPixelType                 m_CostChange;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParallelFastSweepingStep.txx"
#endif

#endif
