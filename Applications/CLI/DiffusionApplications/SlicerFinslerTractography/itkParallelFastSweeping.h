/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParallelFastSweeping.h,v $
  Language:  C++
  Date:      $Date: 2008/02/7 14:28:51 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkParallelFastSweeping_h
#define __itkParallelFastSweeping_h

#include "itkImageToImageFilter.h"
#include "itkParallelFastSweepingStep.h"
#include "itkImage.h"
#include "itkCovariantVector.h"
#include <vector.h>

namespace itk
{
   
/* ******************************************************************************* */
// WORK AROUND BECAUSE NEITHER THE vnl_vector nor the std::vector HAVE A DAMNED
// arg_min() METHOD  (WHATEVER THE DOCUMENTATION SAYS)
class OrderType{
public:
   float        _value;
   unsigned int _position;
   bool operator< ( const OrderType obj ) const
   {
      return( _value < obj._value );
   }
};
/* ******************************************************************************* */
   
/** \class ParallelFastSweeping
 *  \brief Recursively applies iterations of ParallelFastSweepingStep
 *         until convergence based on the total change in the global cost
 *
 */
   
template <class TInputImage, class TOutputImage, class TMaskImage >
class ITK_EXPORT ParallelFastSweeping : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
   /** Convenient typedefs for simplifying declarations. */
   typedef TInputImage                               InputImageType;
   typedef typename InputImageType::Pointer          InputImagePointer;
   typedef typename InputImageType::ConstPointer     InputImageConstPointer;
   typedef TOutputImage                              OutputImageType;
   typedef typename OutputImageType::Pointer         OutputImagePointer;
   typedef typename OutputImageType::ConstPointer    OutputImageConstPointer;
   
   /** Standard class typedefs. */
   typedef ParallelFastSweeping                                 Self;
   typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
   typedef SmartPointer<Self>                                   Pointer;
   typedef SmartPointer<const Self>                             ConstPointer;

   /** Method for creation through the object factory. */
   itkNewMacro(Self);

   /** Run-time type information (and related methods). */
   itkTypeMacro( ParallelFastSweeping, ImageToImageFilter );
  
   /** Image typedef support. */
   typedef typename InputImageType::PixelType   InputPixelType;
   typedef typename OutputImageType::PixelType  OutputPixelType;
   typedef typename OutputImageType::ValueType  ScalarType;
   typedef typename InputImageType::RegionType  InputImageRegionType;
   typedef typename OutputImageType::RegionType OutputImageRegionType;

   typedef typename InputImageType::SizeType    InputSizeType;
   typedef typename InputImageType::IndexType   InputIndexType;
   typedef typename InputImageType::RegionType  InputRegionType;
   typedef typename InputImageType::PointType   InputPointType;
   
   typedef typename OutputImageType::SizeType   OutputSizeType;
   typedef typename OutputImageType::IndexType  OutputIndexType;
   typedef typename OutputImageType::RegionType OutputRegionType;
   typedef typename OutputImageType::PointType  OutputPointType;
   
   /** These types are used to describe the "arrival directions map" (see below) */
   typedef itk::Image<unsigned int,TInputImage::ImageDimension>           ArrivalDirectionsType;
   typedef typename ArrivalDirectionsType::Pointer                        ArrivalDirectionsPointer;
   /** These types are used to store the local costs for each arrival direction */
   typedef itk::VectorImage<OutputPixelType,TInputImage::ImageDimension>  LocalCostType;
   typedef typename LocalCostType::Pointer                                LocalCostPointer;
   /** Type to store the mask for the input, if provided */
   typedef TMaskImage                                                     MaskType;
   typedef typename MaskType::Pointer                                     MaskPointer;
   
   /** Types for the filter performing each fast sweeping step */
   typedef itk::ParallelFastSweepingStep< OutputImageType,
                                          OutputImageType,
                                          InputImageType,
                                          TMaskImage >                    StepFilterType;
   typedef typename StepFilterType::Pointer                               StepFilterPointer;
   
   /** Type to store the arrival directions (neighborhoods) sampled to minimize the local cost */
   typedef itk::CovariantVector<float,TInputImage::ImageDimension>        DirectionType;
   /** Image types for a map of such directions */
   typedef itk::Image< DirectionType, TInputImage::ImageDimension>        DirectionMapType;
   typedef typename DirectionMapType::Pointer                             DirectionMapPointer;
   
   DirectionMapPointer GetOptimumDirectionsMap( void );
   
   //====================================================================================================
   // PARAMETERS TO SET
   void SetMask( MaskPointer mask ){
      m_Mask = mask;
   }
   
   itkSetMacro( MaxIters, unsigned int );
   itkGetMacro( MaxIters, unsigned int );
   
   itkSetMacro( CostFactor, float );
   itkGetMacro( CostFactor, float );
   
   void SetNeighboringDirections( std::vector<DirectionType> dirs )
   {
      m_NeighboringDirections = dirs;
   }
   
   void SetSeedingPoints( std::vector<OutputIndexType> seeds )
   {
      m_SeedingPoints = seeds;
   }
   //====================================================================================================
   // The type of propagation of the solution, used to achieve a certain speed-up
   typedef enum {ALL,CAUSAL,HALF} PropagationType;
   void SetALLPropagationType( void ){
      m_PropagationType = ALL;
   }
   void SetHALFPropagationType( void ){
      m_PropagationType = HALF;
   }
   void SetCAUSALPropagationType( void ){
      m_PropagationType = CAUSAL;
   }
   //====================================================================================================
protected:
   ParallelFastSweeping();
   virtual ~ParallelFastSweeping();
   void PrintSelf(std::ostream& os, Indent indent) const;
   virtual void GenerateInputRequestedRegion();
   void GenerateData();
   //-------------------------------------------------------------
   void CreateMapOfArrivalDirections( void );
   void PrepareNeighborsAndWeights( void );
   void CreateChosenDirectionsPerSwept( void );
   OutputImagePointer CreateCostsMap( void );
   //-------------------------------------------------------------
private:
   ParallelFastSweeping(const Self&); // purposely not implemented
   void operator=(const Self&);       // purposely not implemented
   
   itk::Array2D<unsigned int>   m_Neighbors; // The neighbors to use for each discretized arriving direction
   itk::Array2D<float>          m_Weights;   // The interpolation weights to use for each arriving direction
   /** This is a collection of vectors; each vector of the collection
    corresponds to a sweeping direction for the fast-sweeping algorithm;
    for each sweeping direction, the corresponding vector is a collection
    of indices indicating which sampled spatial directions are actually
    checked to update (or not) the cost at the current voxel. This way we aim
    to exploit the causaity of the solution to the differential equation to
    achieve a net speed-up.*/
   itk::Array<unsigned int>*    m_Chosen;
   /** The following attribute allows to choose the behavior with respect
    to the previous consideration, among the regular behavior (ALL), the
    behavior described above (CAUSAL) or an intermediate behavior where
    we only check neighbors which has been already updated in the current
    swept (HALF).*/
   PropagationType              m_PropagationType;
   /** We keep a scalar image the same size as the input. Each pixel has an 
    unsigned integer value from 0 to the number of directions (i.e. the rows
    of m_Neighbors and m_Weights) minus 1. This value corresponds to the direction
    (or, more generally, the set of neighbors) producing the minimum accumulated
    cost at this voxel. For efficiency purposes, it is the responsibility of the
    user to allocate the required memory for this map and to fix the appropriate
    origin, spacing, regions size... Also, the user must initialize the values
    of the image, that will be iteratively updated by this filter.*/
   ArrivalDirectionsPointer     m_ArrivalDirections;
   /** These are the arrival directions for the fiber bundles. */
   std::vector<DirectionType>   m_NeighboringDirections;
   /** We can use a FA mask to accelerate computations and eliminate
    impossible or unlikely paths:*/
   MaskPointer                  m_Mask;   
   /** A vector where seeding points to perform tractography are stored */
   std::vector<OutputIndexType> m_SeedingPoints;
   
   unsigned int                 m_MaxIters;
   OutputPixelType              m_CostThreshold;
   float                        m_CostFactor;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParallelFastSweeping.txx"
#endif

#endif
