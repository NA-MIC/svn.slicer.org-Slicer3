/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParallelFastSweeping.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkParallelFastSweeping_txx
#define _itkParallelFastSweeping_txx
#include "itkParallelFastSweeping.h"
#include "itkIndexToDirectionImageFilter.h"
#include "itkMatrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_determinant.h"

//===
// DEBUG
#include "itkImageFileWriter.h"
//===

namespace itk
{

/** Constructor */
template <class TInputImage, class TOutputImage, class TMaskImage >
ParallelFastSweeping<TInputImage, TOutputImage, TMaskImage >
::ParallelFastSweeping()
{
   //-----------------------------------------------------------------------
   m_Neighbors.SetSize( 0, TInputImage::ImageDimension );
   m_Weights.SetSize(   0, TInputImage::ImageDimension );
   m_NeighboringDirections.resize(0);
   m_SeedingPoints.resize(0);
   //-----------------------------------------------------------------------
   m_ArrivalDirections = NULL;
   //-----------------------------------------------------------------------
   m_Chosen          = NULL;
   m_PropagationType = HALF;
   //-----------------------------------------------------------------------
   m_MaxIters      = 100;
   m_CostThreshold = itk::NumericTraits<OutputPixelType>::max();
   m_CostFactor    = 0.1;
   m_Mask          = NULL;
   //-----------------------------------------------------------------------
}

/** Destructor */
template <class TInputImage, class TOutputImage, class TMaskImage >
ParallelFastSweeping<TInputImage, TOutputImage, TMaskImage >
::~ParallelFastSweeping()
{
   // Delete c-style allocated memory:
   if( m_Chosen != NULL )
      delete[] m_Chosen;
}
   
/** This filter requires the whole input to work: */
template <class TInputImage, class TOutputImage, class TMaskImage >
typename ParallelFastSweeping<TInputImage, TOutputImage, TMaskImage>::DirectionMapPointer
ParallelFastSweeping<TInputImage, TOutputImage, TMaskImage >
::GetOptimumDirectionsMap( void )
{
   if( !m_ArrivalDirections )
      itkExceptionMacro( << "You must update the filter before calling this method" );
   typedef itk::IndexToDirectionImageFilter<ArrivalDirectionsType,DirectionMapType> MapComputerType;
   typedef typename MapComputerType::Pointer                                        MapComputerPointer;
   MapComputerPointer mapComputer = MapComputerType::New();
   mapComputer->SetNeighboringDirections( m_NeighboringDirections );
   mapComputer->SetInput( m_ArrivalDirections );
   mapComputer->UpdateLargestPossibleRegion();
   return(  mapComputer->GetOutput()  );
}

template< class TInputImage, class TOutputImage, class TMaskImage >
void ParallelFastSweeping< TInputImage, TOutputImage, TMaskImage  >
::CreateMapOfArrivalDirections( void )
{
   m_ArrivalDirections = ArrivalDirectionsType::New();
   m_ArrivalDirections->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
   m_ArrivalDirections->SetOrigin( this->GetInput()->GetOrigin() );
   m_ArrivalDirections->SetSpacing( this->GetInput()->GetSpacing() );
   m_ArrivalDirections->SetDirection( this->GetInput()->GetDirection() );
   m_ArrivalDirections->Allocate();
   m_ArrivalDirections->FillBuffer( m_NeighboringDirections.size() );
}

template< class TInputImage, class TOutputImage, class TMaskImage >
void ParallelFastSweeping< TInputImage, TOutputImage, TMaskImage >
::PrepareNeighborsAndWeights( void )
{
   // Check if the neighboring directions have been set so that we can begin to work:
   if( m_NeighboringDirections.size() == 0 )
      itkExceptionMacro( << "You must fix the set of sampled directions in the sphere" );
   // Check of the map of local costs has been fixed also:
   if( !this->GetInput() )
      itkExceptionMacro( << "You must set the input to this filter, which is a vector image whose voxels are local costs for each sampled direction" );
   // Check that the input has the appropriate number of dimensions:
   if( this->GetInput()->GetVectorLength() != m_NeighboringDirections.size() )
      itkExceptionMacro( << "The number of components of the input has to match the number of neighboring directions" );
   // All seems OK. Let's start computing the neighbors and weights;
   // Allocate memory:
   m_Neighbors.SetSize( m_NeighboringDirections.size(), TInputImage::ImageDimension );
   m_Weights.SetSize(   m_NeighboringDirections.size(), TInputImage::ImageDimension );
   // Compute the number of neighbors in a 3x3x...x3 neighborhood:
   unsigned int numneigh = 1;
   for( unsigned int k=0; k<TInputImage::ImageDimension; ++k ){ numneigh *= 3; }
   // Compute the offsets correspoding to each neighboring position:
   std::vector<DirectionType> offsets(numneigh);
   DirectionType aux;
   aux.Fill( -1 );
   typename std::vector<DirectionType>::iterator it = offsets.begin();
   while( it != offsets.end() ){
      *it = aux;
      for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
         if( aux[d]<1 ){
            aux[d] += 1;
            break;
         }
         else
            aux[d] = -1;
      }
      it++;
   }
   // Temporary vector to allocate distances to each neighbor:
   std::vector<OrderType> distances(numneigh);
   // Proceed:
   for( unsigned int k=0; k<m_NeighboringDirections.size(); ++k ){ // For each sampled direction:
      //----------------------------------------------------------------------------------------------------------------
      // Compute the distance from the corresponding direction to each of the
      // voxels in the 3x3x...x3 neighborhood:
      for( unsigned int n=0; n<numneigh; ++n ){
         distances[n]._value    = ( offsets[n] - m_NeighboringDirections[k] ).GetNorm();
         distances[n]._position = n;
      }
      // Avoid using the central pixel:
      distances[numneigh/2]._value = itk::NumericTraits<float>::max();
      //----------------------------------------------------------------------------------------------------------------
      // Keep the d-closest neighbors:
      std::sort( distances.begin(), distances.end() );
      for( unsigned int d=0; d<TInputImage::ImageDimension; ++d )
         m_Neighbors[k][d] = distances[d]._position;
      //----------------------------------------------------------------------------------------------------------------
      // Compute the interpolation matrix:
      itk::Matrix<float,TInputImage::ImageDimension,TInputImage::ImageDimension> matrix;
      for( unsigned int r=0; r<TInputImage::ImageDimension; ++r ){
         for( unsigned int c=0; c<TInputImage::ImageDimension; ++c )
            matrix[r][c] = offsets[ m_Neighbors[k][c] ][r];
      }
      // To compute the interpolation weights, the three neighbors have to be linearly independent. Make sure of this:
      unsigned int pos = TInputImage::ImageDimension;
      while(   vnl_determinant( matrix.GetVnlMatrix() ) < 1e-3   &&   pos<numneigh  ){
         // Choose the next closest neighbor
         m_Neighbors[k][TInputImage::ImageDimension-1] = distances[pos]._position;
         // Recompute the last column of the matrix:
         for( unsigned int r=0; r<TInputImage::ImageDimension; ++r )
            matrix[r][TInputImage::ImageDimension-1] = offsets[ m_Neighbors[k][TInputImage::ImageDimension-1] ][r];
         // Next neighbor:
         ++pos;
      }
      // Make sure the matrix is not singular:
      if( vnl_determinant( matrix.GetVnlMatrix() ) < 1e-3 )
         itkExceptionMacro( << "Impossible to find suitable weights for some of the neighboring directions" );
      // Invert and compute the weights following J. Melonakos PAMI paper:
      matrix = matrix.GetInverse();
      DirectionType weights = matrix*m_NeighboringDirections[k];
      // Renormalize the weights following J. Melonakos mex code:
      weights.Normalize();
      for( unsigned int d=0; d<TInputImage::ImageDimension; ++d )
         m_Weights[k][d] = ( weights[d]>0 ? weights[d] : -weights[d] ); // Take absolute value of each component
      //----------------------------------------------------------------------------------------------------------------
   }
}

template< class TInputImage, class TOutputImage, class TMaskImage >
typename ParallelFastSweeping< TInputImage, TOutputImage, TMaskImage >::OutputImagePointer
ParallelFastSweeping< TInputImage, TOutputImage, TMaskImage >
::CreateCostsMap( void )
{
   OutputImagePointer costs = OutputImageType::New();
   costs->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
   costs->SetOrigin( this->GetInput()->GetOrigin() );
   costs->SetSpacing( this->GetInput()->GetSpacing() );
   costs->SetDirection( this->GetInput()->GetDirection() );
   costs->Allocate();
   costs->FillBuffer( itk::NumericTraits<OutputPixelType>::max() );
   // Place zeros in seeding points:
   typename std::vector<OutputIndexType>::iterator it;
   it = m_SeedingPoints.begin();
   OutputIndexType idx;
   while( it != m_SeedingPoints.end() ){
      idx = *it;
      bool isInBounds = true;
      for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
         if( idx[d]<0 || idx[d]>=(int)(costs->GetLargestPossibleRegion().GetSize()[d]) )
            isInBounds = false;
      }
      if( isInBounds )
         costs->SetPixel( *it, itk::NumericTraits<OutputPixelType>::Zero );
      it++;
   }
   return costs;
}

/** Initialize the variable m_Chosen: for each sweeping direction, we must
decide which spatial directions among those sampled are actually checked for
cost propagation. We avoid checking those directions violating the causality
of the solution considered at the particular current swept*/
template< class TInputImage, class TOutputImage, class TMaskImage >
void ParallelFastSweeping< TInputImage, TOutputImage, TMaskImage >
::CreateChosenDirectionsPerSwept( void )
{
   // If we have a previous set-up, delete the corresponding memory:
   if( m_Chosen != NULL )
      delete[] m_Chosen;
   // First, determine the number of sweeping directions
   // for this number of dimensions:
   unsigned int numSwepts = 1;
   for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){ numSwepts *= 2; }
   // Allocate the required memory:
   m_Chosen = new itk::Array<unsigned int> [numSwepts];
   // Depending on the preset behavior:
   float signs[TInputImage::ImageDimension];
   switch( m_PropagationType ){
      case ALL:
         // ORIGINAL IMPLEMENTATION: We check each and every spatial
         // directions for all swepts:
         for( unsigned int k=0; k<numSwepts; ++k ){
            // The size of each collection is the total number of sampled
            // spatial directions:
            m_Chosen[k].SetSize( m_NeighboringDirections.size() );
            // Fill this vector with the complete set of indexes:
            for( unsigned int j=0; j<m_NeighboringDirections.size(); ++j )
               m_Chosen[k][j] = j;
         }
         break;
      case CAUSAL:
         // OPTIMIZED IMPLEMENTATION: Only those directions corresponding
         // to previous neighbors according to the causality of the current
         // swept are checked.
         // Create a vector with the signs corresponding to the causality:
         for( unsigned int d=0; d<TInputImage::ImageDimension; ++d )
            signs[d] = -1.0f;
         // For each sweeping direction:
         for( unsigned int k=0; k<numSwepts; ++k ){
            // ------------------------------------------------------------
            // 1- Determine which directions do not violate the causality
            // of the solution in this sweeping direction. They are those
            // whose components have all the same sign as the components
            // in the vector "sign"
            std::vector<unsigned int> chosens;
            chosens.clear();
            for( unsigned int j=0; j<m_NeighboringDirections.size(); ++j ){
               // We choose direction "j" unless one of its signs is wrong:
               bool toBeChosen = true;
               for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
                  if( signs[d]*m_NeighboringDirections[j][d] < -1e-6 )
                     toBeChosen = false;
               }
               // Include the new found direction
               if( toBeChosen )
                  chosens.push_back(j);
            }
            // ------------------------------------------------------------
            // 2- Allocate the required memory for the current itk::Array
            m_Chosen[k].SetSize( chosens.size() );
            // ------------------------------------------------------------
            // 3- Place the indices into the itk::Array
            for( unsigned int j=0; j<chosens.size(); ++j )
               m_Chosen[k][j] = chosens[j];
            // ------------------------------------------------------------
            // 4- Determine the signs for the next sweeping direction
            //    In 2-D, we obtain the four orders: (- -> +), (- -> + )
            //                                       (+ -> -), (- -> + )
            //                                       (- -> +), (+ -> - )
            //                                       (+ -> -), (+ -> - )
            for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
               if( signs[d]<0 ){
                  signs[d] = 1.0f;
                  break;
               }
               else
                  signs[d] = -1.0f;
            }
         }
         break;
      case HALF:
         // Create a vector with the signs corresponding to the causality:
         for( unsigned int d=0; d<TInputImage::ImageDimension; ++d )
            signs[d] = -1.0f;
         // For each sweeping direction:
         for( unsigned int k=0; k<numSwepts; ++k ){
            // ------------------------------------------------------------
            // 1- Determine which directions to consider. They are those whose
            //    last component sign is the same as in the vector of signs;
            //    in case it is zero, we have to check the previous component
            //    and repeat.
            std::vector<unsigned int> chosens;
            chosens.clear();
            for( unsigned int j=0; j<m_NeighboringDirections.size(); ++j ){
               // We choose direction "j" unless one of its signs is wrong:
               bool toBeChosen = true;
               for( int d=TInputImage::ImageDimension-1; d>=0; --d ){
                  if( signs[d]*m_NeighboringDirections[j][d] > 1e-6 )
                     break;
                  else{ // Do we have a wrong sign?
                     if( signs[d]*m_NeighboringDirections[j][d] < -1e-6 ){
                        toBeChosen = false;
                        break;
                     }
                  }
               }
               // Include the new found direction
               if( toBeChosen )
                  chosens.push_back(j);
            }
            // ------------------------------------------------------------
            // 2- Allocate the required memory for the current itk::Array
            m_Chosen[k].SetSize( chosens.size() );
            // ------------------------------------------------------------
            // 3- Place the indices into the itk::Array
            for( unsigned int j=0; j<chosens.size(); ++j )
               m_Chosen[k][j] = chosens[j];
            // ------------------------------------------------------------
            // 4- Determine the signs for the next sweeping direction
            //    In 2-D, we obtain the four orders: (- -> +), (- -> + )
            //                                       (+ -> -), (- -> + )
            //                                       (- -> +), (+ -> - )
            //                                       (+ -> -), (+ -> - )
            for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
               if( signs[d]<0 ){
                  signs[d] = 1.0f;
                  break;
               }
               else
                  signs[d] = -1.0f;
            }
         }
         break;
      default:
         itkExceptionMacro( << "Unsupported sweeping method" );
   }
}
   
template< class TInputImage, class TOutputImage, class TMaskImage >
void ParallelFastSweeping< TInputImage, TOutputImage, TMaskImage >
::GenerateData( void )
{
   //-----------------------------------------------------------------
   // First, create the map of optimal arriving directions (stored in
   // m_ArrivalDirections)
   this->CreateMapOfArrivalDirections();
   //-----------------------------------------------------------------
   // Create auxiliary maps of cots to iterate:
   OutputImagePointer iter0 = this->CreateCostsMap();
   OutputImagePointer iter1 = this->CreateCostsMap();
   //-----------------------------------------------------------------
   // Prepare the neighbors and weights:
   this->PrepareNeighborsAndWeights();
   //-----------------------------------------------------------------
   // Prepare the sets of chosen directions for each swept:
   this->CreateChosenDirectionsPerSwept();
   //-----------------------------------------------------------------
   // Create the filter to perform each step:
   StepFilterPointer step = StepFilterType::New();
   step->SetNeighbors( m_Neighbors );
   step->SetWeights( m_Weights );
   step->SetArrivalDirections( m_ArrivalDirections );
   if( m_Mask )
      step->SetMask( m_Mask );
   // Note the map of local costs is not modified in the subsequent
   // steps, so we can directly pass the input of this filter:
   step->SetLocalCost( this->GetInput() );
   //-----------------------------------------------------------------
   // Iterate until convergence or until a predefined number of iterations:
   bool            mustStop = false;
   unsigned int    iters    = 0;
   OutputPixelType cost0    = itk::NumericTraits<OutputPixelType>::max();
   OutputPixelType cost     = itk::NumericTraits<OutputPixelType>::max();
   while( !mustStop && iters<m_MaxIters ){
      //------------------------------------------------------------------
      // Fix the preferred image direction to split the threads:
      step->SetSplitDirection( iters % TInputImage::ImageDimension );
      //------------------------------------------------------------------
      // Fix the sweeping strategy:
      if( iters==0 ){
         this->SetALLPropagationType();
         this->CreateChosenDirectionsPerSwept();
      }
      if( iters==TInputImage::ImageDimension ){
         this->SetCAUSALPropagationType();
         this->CreateChosenDirectionsPerSwept();
      }
      //else if( iters==TInputImage::ImageDimension )
      //   this->SetHALFPropagationType();
      //else if( iters==2*TInputImage::ImageDimension )
      //   this->SetALLPropagationType();
      step->SetChosen( m_Chosen );
      //------------------------------------------------------------------
      // Fix the input and graft the output to avoid unnecessary memory
      // reallocation:
      if( iters % 2 ){
         iter1->SetRequestedRegion( iter1->GetLargestPossibleRegion() );
         step->SetInput( iter0 );
         step->GraftOutput( iter1 );
         step->Modified();
      }
      else{
         iter0->SetRequestedRegion( iter0->GetLargestPossibleRegion() );
         step->SetInput( iter1 );
         step->GraftOutput( iter0 );
         step->Modified();
      }
      //------------------------------------------------------------------
      // Update and graft the output back
      step->Update();
      if( iters % 2 )
         this->GraftOutput( iter1 );
      else
         this->GraftOutput( iter0 );
      //------------------------------------------------------------------
      // STOP CRITERION:
      // Check if all voxels have a finite value; in this case, we
      // can begin to compute the stop criterion:
      std::cout << "Iteration: " << iters << "; Total change in cost: " << step->GetCostChange() << std::endl;
      /*
      //======================================================================
      //======================================================================
      //======================================================================
      // ONLY FOR DEBUG PURPOSES:
      typedef itk::ImageFileWriter<OutputImageType> WriterType;
      typedef typename WriterType::Pointer          WriterPointer;
      WriterPointer writer = WriterType::New();
      writer->SetInput( step->GetOutput() );
      char name[12];
      sprintf( name, "step%d.mhd", iters );
      writer->SetFileName( name );
      writer->Update();
      //======================================================================
      //======================================================================
      //======================================================================
       */
      if( step->GetCostChange() < m_CostThreshold ){ // Finite cost
         if( cost0 >= m_CostThreshold ){
            // First iteration with finite cost
            cost0 = cost = step->GetCostChange();
         }
         else{
            // There have been previous iterations with finite cost
            cost = step->GetCostChange();
            // Check if the change in the cost is below a fraction
            // of the initial change in the cost:
            mustStop = ( cost <= cost0*m_CostFactor );
         }
      }
      ++iters;
      //------------------------------------------------------------------
   }
   //-----------------------------------------------------------------
   return;
}

/** This filter requires the whole input to work: */
template <class TInputImage, class TOutputImage, class TMaskImage >
void ParallelFastSweeping<TInputImage, TOutputImage, TMaskImage >
::GenerateInputRequestedRegion()
{
   // Call the superclass' implementation of this method
   Superclass::GenerateInputRequestedRegion();
   // Get pointers to the input and output
   InputImagePointer  inputPtr  = const_cast<TInputImage*>( this->GetInput() );
   if ( !inputPtr ){ return; }
   inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
   return;
}

/** Standard "PrintSelf" method */
template <class TInputImage, class TOutput, class TMaskImage >
void ParallelFastSweeping<TInputImage, TOutput, TMaskImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
   Superclass::PrintSelf( os, indent );
   os << indent << "CostThreshold: " << m_CostThreshold << std::endl;
   os << indent << "CostFactor: " << m_CostFactor << std::endl;
   os << indent << "MaxIters: " << m_MaxIters << std::endl;
   
   os << indent << "Neighbors: " << std::endl;
   for( unsigned int r=0; r<m_Neighbors.rows(); ++r ){
      os << indent << indent;
      for( unsigned int c=0; c<m_Neighbors.cols(); ++c )
         os << m_Neighbors[r][c] << indent;
      os << std::endl;
   }
   
   os << indent << "Weights: " << std::endl;
   for( unsigned int r=0; r<m_Neighbors.rows(); ++r ){
      os << indent << indent;
      for( unsigned int c=0; c<m_Neighbors.cols(); ++c )
         os << m_Weights[r][c] << indent;
      os << std::endl;
   }
   
   os << indent << "Neighboring directions: " << std::endl;
   for( unsigned int r=0; r<m_NeighboringDirections.size(); ++r ){
      os << indent << indent << m_NeighboringDirections[r] << std::endl;
   }
   os << std::endl;
   
   os << indent << "Seeding points: " << std::endl;
   for( unsigned int r=0; r<m_SeedingPoints.size(); ++r ){
      os << indent << indent << m_SeedingPoints[r] << std::endl;
   }
   os << std::endl;
   
   if( m_Chosen != NULL ){
      unsigned int numSwepts = 1;
      for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){ numSwepts *= 2; }
      os << indent << "Chosen neighbors for each sweeping direction:" << std::endl;
      for( unsigned int k=0; k<numSwepts; ++k )
         os << indent << indent << "[" << k << "] " << m_Chosen[k] << std::endl;
   }
}

} // end namespace itk

#endif
