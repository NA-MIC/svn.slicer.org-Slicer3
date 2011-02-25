/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParallelFastSweepingStep.txx,v $
  Language:  C++
  Date:      $Date: 2011-01-11 $
  Version:   $Revision: 1.0 $

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkParallelFastSweepingStep_txx
#define __itkParallelFastSweepingStep_txx


#include "itkParallelFastSweepingStep.h"

#include "itkImageRegionIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageDirectionalIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{

template <class TInputImage, class TOutputImage, class TVectorImage, class TMaskImage>
ParallelFastSweepingStep<TInputImage, TOutputImage, TVectorImage, TMaskImage>
::ParallelFastSweepingStep()
{
   //-----------------------------------------------------------------------
   m_CostThreshold = itk::NumericTraits<OutputPixelType>::max();
   //-----------------------------------------------------------------------
   // Precompute the set of offsets for each neighbor:
   unsigned int numneigh = 1;
   for( int i=0; i<TInputImage::ImageDimension; ++i )
      numneigh *= 3;
   m_Offsets.resize( numneigh );
   OutputOffsetType aux;
   aux.Fill( -1 );
   for( unsigned int n=0; n<numneigh; ++n ){
      m_Offsets[n] = aux;
      for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
         if( aux[d]<1 ){
            aux[d] += 1;
            break;
         }
         else
            aux[d] = -1;
      }
   }
   //-----------------------------------------------------------------------
   m_Neighbors.SetSize( numneigh-1, TInputImage::ImageDimension );
   m_Weights.SetSize(   numneigh-1, TInputImage::ImageDimension );
   for( unsigned int k=0; k<numneigh/2; ++k ){
      for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
         m_Weights[k][d]   = 1.0f;
         m_Neighbors[k][d] = k;
      }
   }
   for( unsigned int k=numneigh/2+1; k<numneigh; ++k ){
      for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
         m_Weights[k-1][d]   = 1.0f;
         m_Neighbors[k-1][d] = k;
      }
   }
   //-----------------------------------------------------------------------
   m_Chosen            = NULL;
   //-----------------------------------------------------------------------
   m_ArrivalDirections = NULL;
   m_LocalCost         = NULL;
   m_Mask              = NULL;
   //-----------------------------------------------------------------------
   m_SplitDirection = 1;
   //-----------------------------------------------------------------------
   m_PerThreadCostChange.SetSize( 0 );
   m_CostChange = itk::NumericTraits<OutputPixelType>::max();
}

// Choose a different splitting direction each time (via SetSplitDirection() ) to accelerate
// the propagation of the arrival times and hence improve convergence
template <class TInputImage, class TOutputImage, class TVectorImage, class TMaskImage>
int ParallelFastSweepingStep<TInputImage, TOutputImage, TVectorImage, TMaskImage>
::SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion)
{
   // Get the output pointer
   OutputImageType * outputPtr = this->GetOutput();
   const typename TOutputImage::SizeType& requestedRegionSize = outputPtr->GetRequestedRegion().GetSize();
   
   int splitAxis;
   typename TOutputImage::IndexType splitIndex;
   typename TOutputImage::SizeType splitSize;
   
   // Initialize the splitRegion to the output requested region
   splitRegion = outputPtr->GetRequestedRegion();
   splitIndex = splitRegion.GetIndex();
   splitSize = splitRegion.GetSize();
   
   // The default is to split on the outermost dimension available;
   // instead, we try to split in the direction given by m_SplitDirection
   // and if it is not possible try the default behavior:
   splitAxis = m_SplitDirection;
   if( requestedRegionSize[splitAxis] == 1 ){
      splitAxis = outputPtr->GetImageDimension() - 1;
      while ( requestedRegionSize[splitAxis] == 1 ){
         --splitAxis;
         if (splitAxis < 0){ // cannot split
            itkDebugMacro("  Cannot Split");
            return 1;
         }
      }
   }
   // determine the actual number of pieces that will be generated
   typename TOutputImage::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
   int valuesPerThread = (int)::vcl_ceil(range/(double)num);
   int maxThreadIdUsed = (int)::vcl_ceil(range/(double)valuesPerThread) - 1;
   
   // Split the region
   if (i < maxThreadIdUsed){
      splitIndex[splitAxis] += i*valuesPerThread;
      splitSize[splitAxis] = valuesPerThread;
   }
   if (i == maxThreadIdUsed){
      splitIndex[splitAxis] += i*valuesPerThread;
      // last thread needs to process the "rest" dimension being split
      splitSize[splitAxis] = splitSize[splitAxis] - i*valuesPerThread;
   }
   // set the split region ivars
   splitRegion.SetIndex( splitIndex );
   splitRegion.SetSize( splitSize );
   
   return maxThreadIdUsed + 1;
}

template <class TInputImage, class TOutputImage, class TVectorImage, class TMaskImage>
void ParallelFastSweepingStep<TInputImage, TOutputImage, TVectorImage, TMaskImage>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
   // call the superclass' implementation of this method
   Superclass::GenerateInputRequestedRegion();
   // get pointers to the input and output
   typename Superclass::InputImagePointer inputPtr = const_cast< TInputImage * >( this->GetInput() );
   typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
   
   if ( !inputPtr || !outputPtr )
      return;
   
   // get a copy of the input requested region (should equal the output
   // requested region)
   typename TInputImage::RegionType inputRequestedRegion;
   inputRequestedRegion = inputPtr->GetRequestedRegion();
   
   // pad the input requested region by the operator radius
   InputSizeType radius;
   radius.Fill( 1 );
   inputRequestedRegion.PadByRadius( radius );
   
   // crop the input requested region at the input's largest possible region
   if(   inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() )   ){
      inputPtr->SetRequestedRegion( inputRequestedRegion );
      return;
   }
   else{
      // Couldn't crop the region (requested region is outside the largest
      // possible region).  Throw an exception.
      // store what we tried to request (prior to trying to crop)
      inputPtr->SetRequestedRegion( inputRequestedRegion );
      // build an exception
      InvalidRequestedRegionError e(__FILE__, __LINE__);
      e.SetLocation(ITK_LOCATION);
      e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
      e.SetDataObject(inputPtr);
      throw e;
   }
}

template< class TInputImage, class TOutputImage, class TVectorImage, class TMaskImage>
void ParallelFastSweepingStep< TInputImage, TOutputImage, TVectorImage, TMaskImage>
::BeforeThreadedGenerateData( void )
{
   // Check if we have been told the sampled directions to use at each swept:
   if( !m_Chosen )
      itkExceptionMacro( << "Please, use SetChosen() to provide directions to consider for each swept" );
   // First of all, we copy the input to the output directly, since this filter is
   // recursive and makes use of previous output to compute the current output. We
   // consider the time consumed in this operation is negligible compared to the
   // actual fast-sweeping iterations, and hence do not use multithread:
   typedef itk::ImageRegionConstIterator<TInputImage>  InputIteratorType;
   typedef itk::ImageRegionIterator<TOutputImage>      OutputIteratorType;
   InputIteratorType  iit( this->GetInput(),  this->GetInput()->GetRequestedRegion() );
   OutputIteratorType oit( this->GetOutput(), this->GetInput()->GetRequestedRegion() );
   for( iit.GoToBegin(),oit.GoToBegin(); !oit.IsAtEnd(); ++iit,++oit )
      oit.Set( iit.Get() );
   // Check if the map of local costs is initialized:
   if( !m_LocalCost )
      itkExceptionMacro( << "No map of local costs is available" );
   // Prepare to compute the change in the cost:
   m_PerThreadCostChange.SetSize( this->GetNumberOfThreads() );
   m_PerThreadCostChange.Fill( itk::NumericTraits<OutputPixelType>::Zero );
}

template< class TInputImage, class TOutputImage, class TVectorImage, class TMaskImage>
void ParallelFastSweepingStep< TInputImage, TOutputImage, TVectorImage, TMaskImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId)
{
   // Directional iterator. Note that, contrary to most of the ITK filters,
   // both of them operate on the same image (the output image), hence the
   // need to initialize the output in the BeforeThreadedGenerateData()
   // method directly copying the output.
   typedef itk::ImageDirectionalIteratorWithIndex<OutputImageType>            DirType;
   typedef itk::ImageDirectionalIteratorWithIndex<ArrivalDirectionsType>      DirArrivalType;
   typedef itk::ImageDirectionalConstIteratorWithIndex<InputImageType>        CDirType;
   typedef itk::ImageDirectionalConstIteratorWithIndex<LocalCostType>         LCCDirType;
   typedef itk::ImageDirectionalConstIteratorWithIndex<MaskType>              MaskDirType;
   // Create the iterators:
   DirType     dirw(    this->GetOutput(), outputRegionForThread );
   CDirType    dirr(    this->GetInput(),  outputRegionForThread );
   LCCDirType  dirlc(   m_LocalCost,       outputRegionForThread );
   dirw.SetRadius(0);
   dirr.SetRadius(0);
   dirlc.SetRadius(0);
   dirw.GoToBegin();
   dirr.GoToBegin();
   dirlc.GoToBegin();
   // These two iterators are not always present, depending on:
   //   - if we want to retrieve the map of displacement labels:
   DirArrivalType dira;
   if( m_ArrivalDirections ){
      dira = DirArrivalType( m_ArrivalDirections, outputRegionForThread );
      dira.SetRadius(0);
      dira.GoToBegin();
   }
   //   - if we want to use a mask:
   MaskDirType dirmask;
   if( m_Mask ){
      dirmask = MaskDirType( m_Mask, outputRegionForThread );
      dirmask.SetRadius(0);
      dirmask.GoToBegin();
   }
   //---------------------------------------------------------------------------------------
   // Neighborhood iterators to get the cost in the surrounding voxels:
   typedef itk::ConstantBoundaryCondition<InputImageType>                             ConstantBoundaryConditionR;
   typedef itk::ConstantBoundaryCondition<OutputImageType>                            ConstantBoundaryConditionW;
   typedef itk::ConstNeighborhoodIterator<InputImageType,ConstantBoundaryConditionR>  NeighborRType;
   typedef itk::ConstNeighborhoodIterator<OutputImageType,ConstantBoundaryConditionW> NeighborWType;
   ConstantBoundaryConditionR cbcr;
   ConstantBoundaryConditionW cbcw;
   cbcr.SetConstant( itk::NumericTraits<InputPixelType>::max()  );
   cbcw.SetConstant( itk::NumericTraits<OutputPixelType>::max() );
   OutputSizeType radius;
   radius.Fill(1);
   NeighborRType nItR( radius, this->GetInput(),  outputRegionForThread );
   NeighborWType nItW( radius, this->GetOutput(), outputRegionForThread );
   nItR.OverrideBoundaryCondition( &cbcr );
   nItW.OverrideBoundaryCondition( &cbcw );
   //---------------------------------------------------------------------------------------
   // Determine the buffered region for the output:
   OutputIndexType low  = outputRegionForThread.GetIndex();
   OutputIndexType high = outputRegionForThread.GetIndex() + outputRegionForThread.GetSize();
   // Create buffers to store the pixels for fast access to the neighborhood (note this
   // operation is intensively repeated):
   OutputPixelType* bufferr = new OutputPixelType[ m_Offsets.size() ];
   OutputPixelType* bufferw = new OutputPixelType[ m_Offsets.size() ];
   // For each possible sweeping direction:
   unsigned int sweptDir = 0; // To keep track of the direction inside the loop
   while( !dirw.IsAtEnd() ){
      while( !dirw.IsAtEndOfDirection() ){ // For each voxel in the region
         // If there is a mask and this voxel is not inside it, do nothing
         if( m_Mask ){
            if( !dirmask.Get() ){
               dirw.Set(-1.0f);
               ++dirw;
               ++dirr;
               ++dirlc;
               if( m_ArrivalDirections ){
                  dira.Set( 0 );
                  ++dira;
               }
               ++dirmask;
               continue;
            }
         }
         // Determine the index of this pixel for future use:
         OutputIndexType bound = dirw.GetIndex();
         // Set the position for the neighborhood iterator:
         nItR.SetLocation( bound );
         nItW.SetLocation( bound );
         // Retrieve all the pixels in the neighborhood to the buffer:
         /* THE CODE SHOULD BE RE-WRITTEN IN THE FOLLOWING FASHION:
         typename NeighborRType::NeighborhoodType::Iterator i1 = nItR.GetNeighborhood().Begin();
         typename NeighborWType::NeighborhoodType::Iterator i2 = nItW.GetNeighborhood().Begin();
         unsigned int p = 0;
         while( i1 != nItR.GetNeighborhood().End() ){
            if( p==m_Offsets.size() )
               std::cerr << bound << std::endl;
            bufferr[p] = *i1;
            bufferw[p] = *i2;
            ++p; ++i1; ++i2;
         }
         BUT SUCH CODE CRASHES FOR SOME REASON
         */
         for( unsigned int p=0; p<m_Offsets.size(); ++p ){
            bufferr[p] = nItR.GetPixel(p);
            bufferw[p] = nItW.GetPixel(p);
         }
         /** --------------------------------------------------------------------------------------- */
         /** --------------------------------------------------------------------------------------- */
         /** --------------------------------------------------------------------------------------- */
         /*
         // FORMER IMPLEMENTATION WITHOUT A MASK:
         // If this pixel is infinite, the cost has to be fixed also to infinite change, so that
         // additional iterations are required
         if( dirw.Get() >= m_CostThreshold )
            m_PerThreadCostChange[threadId] = itk::NumericTraits<OutputPixelType>::max();
         */
         // BUT IF WE CONSIDER WE MAY HAVE A MASK, WE CANNOT DO SUCH THING; IN CASE A WHITE ISLE
         // SURROUNDED BY BLACK VOXELS (OR THE IMAGE BOUNDARY) WHERE PRESENT, A FINITE COST WOULD
         // NEVER PROPAGATE TO THE VOXELS IN THE AISLE, FORCING THE ALGORITHM TO EXECUTE UNTIL THE
         // MAXIMUM NUMBER OF ITERATIONS IS REACHED. INSTEAD, WE ONLY SET THE PER THREAD COST TO
         // INFINITE VALUE IF THE ALUE OF THE VOXEL IS ACTUALLY UPDATED. AT EACH ITERATION, AT LEAST
         // ONE INFINITE VOXEL IN THE IMAGE THAT HAS TO CHANGE WILL BE UPDATED, HENCE THE PROCEDURE
         // IS CORRECT   
         /** --------------------------------------------------------------------------------------- */
         /** --------------------------------------------------------------------------------------- */
         /** --------------------------------------------------------------------------------------- */
         // We have to check all possible sampling directions of the sphere. For each of them:
         //   1- Decide which neighbors in the 1x1x...x1 neigborhood we use. The number of neighbors
         //      is the same as the image dimension. They can be precomputed for each sampled direction,
         //      and it is the responsibility of the calling filter to set them via the SetNeighbors()
         //      method (for efficiency reasons).
         //   2- Interpolate the arrival times for these neighbors, with weights that can be also
         //      precomputed and must be set by the calling filter via SetWeights().
         //   3- From the interpolated arrival times and the local cost for the corresponding direction,
         //      compute the local arrival time for this particular direction. Update the arrival time
         //      only if the result is smaller than the previous one.
         OutputPixelType minimumCost = itk::NumericTraits<OutputPixelType>::max();
         unsigned int    arg         = m_Neighbors.rows();
         // For each discretized direction among those to be checked for this particular
         // sweeping direction:
         //for( unsigned int g=0; g<m_Neighbors.rows(); ++g ){
         for( unsigned int g_s=0; g_s<m_Chosen[sweptDir].Size(); ++g_s ){
            // This is the actual direction in the global list (note all arrays are
            // indexed with respect to g, which is the global one)
            unsigned int g = m_Chosen[sweptDir][g_s];
            // Interpolate the arrival times from the neighbors corresponding to the
            // current arrival direction:
            float cost = itk::NumericTraits<float>::Zero;   // Interpolated value
            float norm = itk::NumericTraits<float>::Zero;   // Normalization for the interpolation
            for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){ // For each neighbor
               // 1- Check if this neighbor is inside the region for thread. Otherwise,
               //    we have to pick it from the input (not the output) to keep the algorithm
               //    thread-safe
               OutputPixelType val = itk::NumericTraits<OutputPixelType>::max();
               OutputIndexType pos = bound + m_Offsets[ m_Neighbors[g][d] ];
               bool isInside = true;
               for( unsigned int k=0; k<TInputImage::ImageDimension; ++k ){
                  if( pos[k]<low[k] || pos[k]>=high[k] )
                     isInside = false;
               }
               if( isInside )
                  val = bufferw[ m_Neighbors[g][d] ];
               else
                  val = bufferr[ m_Neighbors[g][d] ];
               // 2- Check if the cost is finite:
               bool isInf = ( val >= m_CostThreshold );
               // 3- If the neighbor is finite, we can average
               if( !isInf && val>=0 ){
                  cost += ( m_Weights[g][d] * (float)val );
                  norm += m_Weights[g][d];
               }
            }
            // Compute the arrival time from the local cost
            if( norm>1e-5 ){ // I.e., if any of the neighbors has a finite arrival cost
               cost = ( cost + dirlc.Get()[g] )/norm;
               // Keep the minimum cost and the corresponding argument
               if( cost<minimumCost ){
                  minimumCost = static_cast<OutputPixelType>( cost );
                  arg         = g;
               }
            }
         }
         // In minimumCost we have the minimum arrival time corresponding to the optimum
         // direction from which the current voxel is reached. The overall cost has to be
         // updated only if this arrival time is smaller than the original one:
         if( minimumCost < dirw.Get() ){
            /*
            // FORMER IMPLEMENTATION WITHOUT FA MASK:
            // Update the cost per thread, if its not already infinite:
            if( m_PerThreadCostChange[threadId] < m_CostThreshold )
               m_PerThreadCostChange[threadId] += ( dirw.Get() - minimumCost );
            */
            /** ---------------------------------------------------------------------------------- */
            // NEW IMPLEMENTATION TO USE WITH A FA MASK:
            if( dirw.Get() >= m_CostThreshold )
               m_PerThreadCostChange[threadId] = itk::NumericTraits<OutputPixelType>::max();
            else if( m_PerThreadCostChange[threadId] < m_CostThreshold )
               m_PerThreadCostChange[threadId] += ( dirw.Get() - minimumCost );
            /** ---------------------------------------------------------------------------------- */
            dirw.Set( minimumCost );
            // If the map of arrival directions is present, it must be updated too
            if( m_ArrivalDirections )
               dira.Set( arg );
         }
         ++dirw;
         ++dirr;
         ++dirlc;
         if( m_ArrivalDirections )
            ++dira;
         if( m_Mask )
            ++dirmask;
      }
      ++sweptDir;
      dirw.NextDirection();
      dirr.NextDirection();
      dirlc.NextDirection();
      if( m_ArrivalDirections )
         dira.NextDirection();
      if( m_Mask )
         dirmask.NextDirection();
   }
   // Delete the allocated buffer:
   delete[] bufferr;
   delete[] bufferw;
}
   
template< class TInputImage, class TOutputImage, class TVectorImage, class TMaskImage>
void ParallelFastSweepingStep< TInputImage, TOutputImage, TVectorImage, TMaskImage>
::AfterThreadedGenerateData( void )
{
   // Compute the global change in the cost
   m_CostChange = itk::NumericTraits<OutputPixelType>::Zero;
   for( int k=0; k<this->GetNumberOfThreads(); ++k ){
      if( m_PerThreadCostChange[k] >= m_CostThreshold ){
         m_CostChange = itk::NumericTraits<OutputPixelType>::max();
         break;
      }
      m_CostChange += m_PerThreadCostChange[k];
   }
}
   
template <class TInputImage, class TOutput, class TVectorImage, class TMaskImage>
void ParallelFastSweepingStep<TInputImage, TOutput, TVectorImage, TMaskImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
   Superclass::PrintSelf( os, indent );
   os << indent << "Split Direction: " << m_SplitDirection << std::endl;
   os << indent << "Cost change: "     << m_CostChange     << std::endl;
   os << indent << "Cost threshold: "  << m_CostThreshold  << std::endl;
   os << indent << "Neighbors:      "  << std::endl << indent << indent << m_Neighbors << std::endl;
   os << indent << "Weights:        "  << std::endl << indent << indent << m_Weights   << std::endl;
   os << indent << "Offsets:        "  << std::endl;
   for( unsigned int i=0; i<m_Offsets.size(); ++i )
      os << indent << indent << m_Offsets[i] << std::endl;
}


} // end namespace itk


#endif
