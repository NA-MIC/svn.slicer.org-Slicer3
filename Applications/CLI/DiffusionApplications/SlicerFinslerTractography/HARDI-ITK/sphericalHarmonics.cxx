#ifndef _sphericalHarmonics_cxx
#define _sphericalHarmonics_cxx

#include "sphericalHarmonics.h"
#include "math.h"

namespace shmaths
{

/** Cumulative product. It returns factorial(end)/factorial(start) or, alternatively,
 (start+1)·(start+2)·...·(end-1)·end. Note that necessarily end>start; otherwise, it
 returns simply 1. This function can be used only with SMALL unsigned integer numbers
 since it uses integer arithmetic. Otherwise, it could produce overflow.
 
 Last version: 07/05/2010
 */
unsigned long cumprod( const unsigned int start, const unsigned int end )
{
   if( start>=end )
      return 1;
   else
      return( cumprod(start,end-1) * end );
}

/** Compute the value of the associated Legendre polynomials of degrees l = 0, 1, ... L
 at the desired point x. Such point must satisfy the condition -1 <= x <= 1, and NO CHECKING
 IS PERFORMED TO THIS RESPECT. The values computed are stored in "buffer" which sould
    be initialized with a call of the form:
    
    float* buffer = allocateBufferForAssociatedLegendrePolynomials( L );
 
 for the same value of L as used in this function. This buffer has to be erased later
 on by means of the delete[] operator. Note that for each degree l (2l+1) polynomials,
 corresponding to orders m = -l, ..., 0, ... l have to be evaluated. Hence, the total
 amount of polynomials to be evaluated grow to (L+1)^2.
 
 To compute the desired values, we used the following recursive formula (m>1):
 
 sqrt(1-x^2) P_l^{m+1}(x) = (l-m)·x·P_l^m(x) - (l+m)·P_(l-1)^m(x)
 
 and for m=0 we have that P_l^0(x) = P_l(x), the standard Legendre polynomial. Thsu, we 
 can use the function computeLegendrePolynomials to start the recursive formula.
 
 For negative values of m, we can use the relation:
 
 P_l^{-m}(x) = (-1)^m·(l-m)!/(l+m)!·P_l^m(x)
 
 The function cumprod implemented above can be used for this computation.
 
 Finally, note that a singularity may arise for x=1 in the recursive formula. It is easy
 to check that in fact teh value of all associated Legendre polynomials with l,m>0 is 0
 for x=+-1, so we use this property to crop the values of the polynomials
 
 Last version: 07/05/2010
*/
void computeAssociatedLegendrePolynomials( const double x, const unsigned int L, double* buffer )
{
   // First, compute the values of the standard Legendre polynomials, which are to be used
   // to initalized the recursive formula (they are the values for m=0).
   double* aux = allocateBufferForLegendrePolynomials( L );
   computeLegendrePolynomials( x, L, aux );
   // Compute the associated polynomials recursively
   double sqx = x*x;
   if( (1-sqx)>1e-10 ){ // There are no numerical issues, so we can compute the value normally
      // For degree 0, we have the first position and only one polynomial; note that L
      // is unsigned int, hence this position is always defined and no bound checking is needed:
      buffer[0] = aux[0];
      for( unsigned int l=1; l<=L; ++l ){
         // The position of degree l, m=0 in the actual indexing of the  buffer (current l index):
         unsigned int cli = l*l + l;
         // The position of degree l-1, m=0 in the actual indexing of the buffer (previous l index):
         unsigned int pli = (l-1)*(l-1) + (l-1);
         // The value for degree l and m=0 is the same as the value of the standard Legendre
         // polynomial, which has been already computed:
         buffer[cli] = aux[l];
         // Now, use the recursive rule to compute the values for m=1..l
         for( unsigned int m=1; m<=l; ++m )
            buffer[cli+m] = (   ((int)l-(int)m+1)*x*buffer[(int)(cli+m)-1]   -   ((int)(l+m)-1)*buffer[(int)(pli+m)-1]   ) / ::sqrt(1-sqx);
         // For negative values of m, we can perform the computation based on the curresponding
         // positive values of m:
         int sign = -1; // This is an auxiliar variable to compute (-1)^m
         for( unsigned int m=1; m<=l; ++m, sign*=(-1) )
            buffer[(int)cli-(int)m] = (sign*buffer[cli+m])/cumprod( (int)l-(int)m, l+m );
         // We're done!!!
      }
   }
   else{ // x is too close to +- 1.0, so we have to crop the value of the polynomial:
      unsigned int pos = 0; // Auxiliar variable
      for( unsigned int l=0; l<=L; ++l ){ // For each order:
         // Fill negative values of m:
         for( int m=-(int)l; m<0; ++m )
            buffer[pos++] = 0.0f;
         // Fix the value for m=0
         buffer[pos++] = aux[l];
         // Fill positive values of m:
         for( int m=1; m<=(int)l; ++m )
            buffer[pos++] = 0.0f;
      }
   }
   // Delete the auxiliar buffer where standard Legendre polynomials are stored:
   delete[] aux;
   return;
}

/** Compute the value of the Legendre polynomials of degrees l = 0, 1, ... L at the 
 desired point x. Such point must satisfy the condition -1 <= x <= 1, and NO CHECKING
 IS PERFORMED TO THIS RESPECT. The values computed are stored in "buffer" which sould
 be initialized with a call of the form:
 
 float* buffer = allocateBufferForLegendrePolynomials( L );
 
 for the same value of L as used in this function. This buffer has to be erased later 
 on by means of the delete[] operator.
 
 To compute the desired values, we used the following recursive formula:
 
   P_0(x)  = 1;
   P_1(x)  = x;
   lP_l(x) = (2l-1)xP_(l-1)(x) - (l-1)P_(l-2)(x)
 
 Last version: 07/05/2010
*/
void computeLegendrePolynomials( const double x, const unsigned int L, double* buffer )
{
   for( unsigned int l=0; l<=L; ++l ){
      if( l==0 )
         buffer[l] = 1;
      else if( l==1 )
         buffer[l] = x;
      else
         buffer[l] = ( (2*l-1)*x*buffer[l-1] - (l-1)*buffer[l-2] ) / l;
   }
}

/** Allocate the buffer necessary to store the evaluations of all associated
 Legendre polynomials of order l = 0, 1, ..., L at a given real number -1 <= x <= 1.
 For each degree l, possible orders range from m=-l to m=l, so we have to allocate
 $\sum_{l=0}^{L}(2*l+1) = (L+1)^2$ doubles.
 
 THE MEMORY ALLOCATED BY THIS FUNCTION MUST BE DELETED EXTERNALLY WITH THE
 delete[] OPERATOR
 
 Last version: 07/05/2010
 */
double* allocateBufferForAssociatedLegendrePolynomials( const unsigned int L )
{
   return (new double[(L+1)*(L+1)]);
}

/** Allocate the buffer necessary to store the evaluations of all Legendre 
 polynomials of order l = 0, 1, ..., L at a given real number -1 <= x <= 1.
 We have to allocate $\sum_{l=0}^{L}1 = (L+1)$ doubles.
 
 THE MEMORY ALLOCATE WITH THIS FUNCTION MUST BE DELETED EXTERNALLY WITH THE
 delete[] OPERATOR
 
 Last version: 07/05/2010
 */
double* allocateBufferForLegendrePolynomials( const unsigned int L )
{
   return (new double[L+1]);
}

   
/** Compute the value of the associated Legendre polynomials of even degrees l = 0, 2, ... L
 at the desired point x. Such point must satisfy the condition -1 <= x <= 1, and NO CHECKING
 IS PERFORMED TO THIS RESPECT. The values computed are stored in "buffer" which sould
 be initialized with a call of the form:
 
 float* buffer = allocateBufferForEvenAssociatedLegendrePolynomials( L );
 
 for the same value of L as used in this function. This buffer has to be erased later
 on by means of the delete[] operator. Note that for each degree l (2l+1) polynomials,
 corresponding to orders m = -l, ..., 0, ... l have to be evaluated. Hence, the total
 amount of polynomials to be evaluated grow to (L+1)·(L/2+1) (only even degrees are
 stored). L IS ASSUMED TO BE EVEN, AND NO CHECKING IS PERFORMED TO THIS RESPECT.
 
 This function is based on computeAssociatedLegendrePolynomials; in fact, the values for all
 degrees (both even and odd) are computed, and only those corresponding to even orders are
 stored in the final buffer. This is because we need a recursive formula to compute the value
 of the polynomials, involving even and odd degree polyomials. Note that this function is
 provided only for convenience, since it is very useful in problems showing spherical symmetry.
 
 Last version: 07/05/2010
 */
void computeEvenAssociatedLegendrePolynomials( const double x, const unsigned int L, double* buffer )
{
   // Compute all the associated Legendre polynomials, both for even and odd orders:
   double* aux = allocateBufferForAssociatedLegendrePolynomials( L );
   computeAssociatedLegendrePolynomials( x, L, aux );
   // We have to keep only the polynomials associated to even orders:
   for( unsigned int l=0; l<=L/2; ++l ){
      // The position of degree 2l, m=0 in the actual indexing of the auxiliar buffer (auxiliar l index):
      unsigned int ali = (2*l)*(2*l) + (2*l);
      // The position of degree 2l, m=0 in the actual indexing of the final buffer (final buffer l index):
      unsigned int fli = 2*l*l + l;
      // Now, we can copy from the auxiliar buffer to the final buffer for each order m:
      for( int m=-(int)(2*l); m<=(int)(2*l); ++m )
         buffer[(int)fli+m] = aux[(int)ali+m];
   }
   // Delete the buffer with all the associated polynomials:
   delete[] aux;
}

/** Compute the value of the Legendre polynomials of eben degrees l = 0, 2, 4, ... L
 at the desired point x. Such point must satisfy the condition -1 <= x <= 1, and NO CHECKING
 IS PERFORMED TO THIS RESPECT. The values computed are stored in "buffer" which sould
 be initialized with a call of the form:
 
 float* buffer = allocateBufferForEvenLegendrePolynomials( L );
 
 for the same value of L as used in this function. This buffer has to be erased later 
 on by means of the delete[] operator. L IS ASSUMED TO BE EVEN, AND NO CHECKING IS PERFORMED
 TO THIS RESPECT.
 
 This function is based on computeLegendrePolynomials; in fact, the values for all degrees
 (both even and odd) are computed, and only those corresponding to even orders are stored
 in the final buffer. This is because we need a recursive formula to compute the value of the 
 polynomials, involving even and odd degree polyomials. Note that this function is provided
 only for convenience, since it is very useful in problems showing spherical symmetry.
 
 Last version: 07/05/2010
 */
   
void computeEvenLegendrePolynomials( const double x, const unsigned int L, double* buffer )
{
   // Use the general purpose function to compute the values for all orders:
   double* aux = allocateBufferForLegendrePolynomials( L );
   computeLegendrePolynomials( x, L, aux );
   // Store the values only for even degrees:
   for( unsigned int l=0; l<=L/2; ++l )
      buffer[l] = aux[2*l];
   // Delete auxiliar buffer:
   delete[] aux;
}

/** Allocate the buffer necessary to store the evaluations of all even degree associated
 Legendre polynomials of order l = 0, 2, ..., L at a given real number -1 <= x <= 1.
 For each degree l, possible orders range from m=-l to m=l, so we have to allocate
 $\sum_{l=0}^{L/2}(2*(2l)+1) = (L+1)·(L/2+1)$ doubles.
 
 L IS ASSUMED TO BE ODD, AND NO OTHER CHECKING IS DONE
 
 THE MEMORY ALLOCATE WITH THIS FUNCTION MUST BE DELETED EXTERNALLY WITH THE
 delete[] OPERATOR
 
 Last version: 07/05/2010
 */
double* allocateBufferForEvenAssociatedLegendrePolynomials( const unsigned int L )
{
   return (new double[(L+1)*(L/2+1)]);
}

/** Allocate the buffer necessary to store the evaluations of all Legendre 
 polynomials of even orders l = 0, 2, ..., L at a given real number -1 <= x <= 1.
 We have to allocate $\sum_{l=0}^{L/2}1 = (L/2+1)$ doubles.
 
 L IS ASSUMED TO BE ODD, AND NO OTHER CHECKING IS DONE
 
 THE MEMORY ALLOCATE WITH THIS FUNCTION MUST BE DELETED EXTERNALLY WITH THE
 delete[] OPERATOR
 
 Last version: 07/05/2010
 */
double* allocateBufferForEvenLegendrePolynomials( const unsigned int L )
{
   return (new double[L/2+1]);
}
   
/** Return the number of Legendre polynomials up to degree L */
unsigned int getNumberOfLegendrePolynomials( const unsigned int L )
{
   return (L+1);
}

/** Return the number of associated Legendre polynomials up degree order L */
unsigned int getNumberOfAssociatedLegendrePolynomials( const unsigned int L )
{
   return ( (L+1)*(L+1) );
}

/** Return the number of Legendre polynomials up to degree L (only even degrees).
 L IS ASSUMED TO BE EVEN, AND NO CHECKING IS PERFORMED*/
unsigned int getNumberOfEvenLegendrePolynomials( const unsigned int L )
{
   return (L/2+1);
}

/** Return the number of associated Legendre polynomials up to degree L (only even degrees).
 L IS ASSUMED TO BE EVEN, AND NO CHECKING IS PERFORMED*/
unsigned int getNumberOfEvenAssociatedLegendrePolynomials( const unsigned int L )
{
   return ( (L+1)*(L/2+1) );
}
   
/** Compute the value of Legendre polynomial of degree l at the point x, which
 is assumed to lay in the range [-1,1] (NO CHECKING IS PERFORMED). Note that this
 function internally uses computeLegendrePolynomials, so the values for all degrees
 r<l have to be computed each time the function is called. Consider using instead:
 
 double computeP_l( const unsigned int l, const double* buffer );
 
 where buffer has to be computed with the aforementioned computeLegendrePolynomials
 for the same value of x. For succesive calls with the same argument x, precomputing
 the buffer can considerably accelerate the execution.
*/
double computeP_l( const double x, const unsigned int l )
{
   double*  buffer = allocateBufferForLegendrePolynomials( l );
   computeLegendrePolynomials( x, l, buffer );
   double   value  = computeP_l( l, buffer );
   delete[] buffer;
   return (value);
}

/** Compute the value of associated Legendre polynomial of degree l ans order m
 at the point x, which is assumed to lay in the range [-1,1] (NO CHECKING IS
 PERFORMED). Note that this function internally uses
 computeAssociatedLegendrePolynomials, so the values for all degrees r<l have
 to be computed each time the function is called. Consider using instead:
 
 double computeP_l( const unsigned int l, const int m, double* buffer );

 where buffer has to be computed with the aforementioned
 computeAssociatedLegendrePolynomials for the same value of x. For succesive
 calls with the same argument x, precomputing the buffer can considerably accelerate
 the execution.
 
 Note that m is assumed to lay in the range [-l,l], and no bound checking is
 performed.
 */
double computeP_l_m( const double x, const unsigned int l, const int m )
{
   double*  buffer = allocateBufferForAssociatedLegendrePolynomials( l );
   computeAssociatedLegendrePolynomials( x, l, buffer );
   double   value = computeP_l_m( l, m, buffer );
   delete[] buffer;
   return (value);
}

/** Compute the value of Legendre polynomial of degree l at the point x, which
 is assumed to lay in the range [-1,1] (NO CHECKING IS PERFORMED). The buffer
 argument can be obtained via the computeLegendrePolynomials function.
 THE USER IS RESPONSIBLE OF EXTERNALLY ALLOCATING, COMPUTING AND DELETING
 THE BUFFER
*/
double computeP_l( const unsigned int l, const double* buffer )
{
   return (buffer[l]);
}

/** Compute the value of associated Legendre polynomial of degree l ans order m
 at the point x, which is assumed to lay in the range [-1,1] (NO CHECKING IS
 PERFORMED). The buffer argument can be obtained via the
 computeAssociatedLegendrePolynomials function. THE USER IS RESPONSIBLE OF
 EXTERNALLY ALLOCATING, COMPUTING AND DELETING THE BUFFER
 
 Note that m is assumed to lay in the range [-l,l], and no bound checking
 is performed.
 */
double computeP_l_m( const unsigned int l, const int m, const double* buffer )
{
   return (   buffer[ (int)(l*(l+1)) + m ]   );
}

/** Compute the real-valued Spherical Harmonics basis function of degree l and
 order m for a point given, in spherical coordinates, as:
 
   (r=1,theta,phi)
 
 Note that physics convention is used, so that the unit vectors on the
 canonycal cartesian basis are given by:
 
  (1,0,0) -> ( r=1, theta=pi/2, phi=0         )
  (0,1,0) -> ( r=1, theta=pi/2, phi=pi/2      )
  (0,0,1) -> ( r=1, theta=0,    phi=undefined )
 
 This is a modified, real-valued basis; instead of complex exponentials,
 cosines (for m>0) and sines (for m<0) are used.
 */
double computeY_l_m( const double phi, const double theta, const unsigned int l, const int m )
{
   if( m==0 )
      return (   sqrt( (2*l+1)/(4*PI) ) * computeP_l( cos(theta), l )   );
   else if( m>0 )
      return (   sqrt( (2*l+1)/(2*PI)/cumprod((int)l-m,(int)l+m) ) * computeP_l_m( cos(theta), l, m )   ) * sin( m*phi );
   else
      return (   sqrt( (2*l+1)/(2*PI)*cumprod((int)l+m,(int)l-m) ) * computeP_l_m( cos(theta), l, m )   ) * cos( m*phi );      
}

/** Compute the whole matrix of spherical harmonics for least squeares fitting.
 The resulting matrix has size NxH. Hence, N represents the number of points over
 the unit-radius sphere for which the real SH-basis is evaluated (therefore, both
 theta and phi have to be length-N vectors). Each row of the resulting matrix sh
 represents the evaluation of the first H basis functions for the point defined
 by spherical coordinates (r=1, theta, phi) with physics convention. The first
 column represents (l=0,m=0). The second one is (l=1,m=-1); then (l=1,m=0), 
 (l=1,m=1),(l=2,m=-2), and so on (spherical harmonics are provided in the same
 order as the associated Legendre polynomials).
 */
void computeSHMatrix( const unsigned int N, const double* theta, const double* phi, const unsigned int L, SHMatrixType& sh )
{
   // Set the proper size for the resulting matrix:
   sh.set_size( N, getNumberOfAssociatedLegendrePolynomials(L) );
   // Allocate the buffer to precompute the associanted Legendre polynomials:
   double* buffer = allocateBufferForAssociatedLegendrePolynomials( L );
   // Now, for each point in the unit sphere (and hence, for each row of the matrix):
   for( unsigned int n=0; n<N; ++n ){
      // Get the point to evaluate the polynomials:
      double x0 = cos( theta[n] );
      // Precompute the buffer of associated Legendre polynomials:
      computeAssociatedLegendrePolynomials( x0, L, buffer );
      // Now, compute the SH values and store them in the same buffer:
      unsigned int pos=0; // Auxiliar position in the buffer
      for( unsigned int l=0; l<=L; ++l ){  // For each degree
         for( int m=-(int)l; m<0; ++m ){  // For each negative order
            buffer[pos] = sqrt( (2*l+1)/(2*PI)*cumprod((int)l+m,(int)l-m) ) * buffer[pos] * cos( m * phi[n] );
            ++pos;
         }
         // For m=0:
         buffer[pos] = sqrt( (2*l+1)/(4*PI) ) * buffer[pos];
         ++pos;
         for( int m=1; m<=(int)l; ++m ){  // For each positive order
            buffer[pos] = sqrt( (2*l+1)/(2*PI)/cumprod((int)l-m,(int)l+m) ) * buffer[pos] * sin( m * phi[n] );
            ++pos;
         }
      }
      // Now the buffer contains the whole n-th row of the resulting matrix. We can insert it directly:
      sh.set_row( n, buffer );
   }
   // Delete the buffer:
   delete[] buffer;
}

   
/** Compute the whole matrix of spherical harmonics for least squeares fitting.
 The resulting matrix has size NxH. Hence, N represents the number of points over
 the unit-radius sphere for which the real SH-basis is evaluated (therefore, both
 theta and phi have to be length-N vectors). Each row of the resulting matrix sh
 represents the evaluation of the first H basis functions for the point defined
 by spherical coordinates (r=1, theta, phi) with physics convention. The first
 column represents (l=0,m=0). The second one is (l=1,m=-1); then (l=1,m=0), 
 (l=1,m=1),(l=2,m=-2), and so on (spherical harmonics are provided in the same
 order as the associated Legendre polynomials).
 
 As opposed to computeSHMatrix, this function consideres only even degrees of
 the spherical harmonics. This implies that only functions with radial symmetry
 can be represented.
 */
void computeSHMatrixSymmetric( const unsigned int N, const double* theta, const double* phi, const unsigned int L, SHMatrixType& sh )
{
   // Set the proper size for the resulting matrix:
   sh.set_size( N, getNumberOfEvenAssociatedLegendrePolynomials(L) );
   // Allocate the buffer to precompute the associanted Legendre polynomials:
   double* buffer = allocateBufferForEvenAssociatedLegendrePolynomials( L );
   // Now, for each point in the unit sphere (and hence, for each row of the matrix):
   for( unsigned int n=0; n<N; ++n ){
      // Get the point to evaluate the polynomials:
      double x0 = cos( theta[n] );
      // Precompute the buffer of associated Legendre polynomials:
      computeEvenAssociatedLegendrePolynomials( x0, L, buffer );
      // Now, compute the SH values and store them in the same buffer:
      unsigned int pos=0; // Auxiliar position in the buffer
      for( unsigned int l=0; l<=L/2; ++l ){  // For each even degree
         for( int m=-(int)(2*l); m<0; ++m ){  // For each negative order
            buffer[pos] = sqrt( (4*l+1)/(2*PI)*cumprod((int)(2*l)+m,(int)(2*l)-m) ) * buffer[pos] * cos( m * phi[n] );
            ++pos;
         }
         // For m=0:
         buffer[pos] = sqrt( (4*l+1)/(4*PI) ) * buffer[pos];
         ++pos;
         for( int m=1; m<=(int)(2*l); ++m ){  // For each positive order
            buffer[pos] = sqrt( (4*l+1)/(2*PI)/cumprod((int)(2*l)-m,(int)(2*l)+m) ) * buffer[pos] * sin( m * phi[n] );
            ++pos;
         }
      }
      // Now the buffer contains the whole n-th row of the resulting matrix. We can insert it directly:
      sh.set_row( n, buffer );
   }
   // Delete the buffer:
   delete[] buffer;
}

   
/** This function generates a diagonal matrix (eig) whose diagonal entries are
 the eigenvalues of the Spherical Harmonics basis functions associated to the
 Laplace-Beltrami operator (i.e., the part of the Laplacian depending only on
 the angular coordinates (theta,phi). The eigenvalues for the SH functions up
 to degree L are computed and stored in the diagonal of eig with the usual 
 order: (l=0,m=0), (l=1,m=-1), ..., (l=1,m=1), (l=2,m=-2), ..., (l=2,m=2), ...
 (l=L,m=-L), ...(l=L,m=L).
 
 Interestingly, the associated eigenvalue for each l is the same independently
 on the value of m, and it is computed as: \lambda_{l,m} = -l*(l+1).
 */
void computeSHEigMatrix( const unsigned int L, SHMatrixType& eig )
{
   // Set the proper size:
   eig.set_size( getNumberOfAssociatedLegendrePolynomials(L), getNumberOfAssociatedLegendrePolynomials(L) );
   // Fill with zeros:
   eig.fill( 0.0f );
   // Set the proper values on the diagonal of the matrix:
   unsigned int pos = 0; // Auxiliar position counter
   for( unsigned int l=0; l<=L; ++l ){ // For each degree
      for( int m=-(int)l; m<=(int)l; ++m, ++pos )
         eig(pos,pos) = -(double)( l*(l+1) );
   }
}
   
/** This function generates a diagonal matrix (eig) whose diagonal entries are
 the eigenvalues of the Spherical Harmonics basis functions associated to the
 Laplace-Beltrami operator (i.e., the part of the Laplacian depending only on
 the angular coordinates (theta,phi). The eigenvalues for the SH functions up
 to degree L are computed and stored in the diagonal of eig with the usual
 order: (l=0,m=0), (l=1,m=-1), ..., (l=1,m=1), (l=2,m=-2), ..., (l=2,m=2), ...
 (l=L,m=-L), ...(l=L,m=L).
 
 Interestingly, the associated eigenvalue for each l is the same independently
 on the value of m, and it is computed as: \lambda_{l,m} = -l*(l+1)
 
 As opposed to computeSHEigMatrix, this function is restricted to even degrees
 of the SH basis, and hence it is only usuful for fucntions defined over the
 unit sphere showing radial symmetry. L IS ASSUME DTO BE EVEN AND NO CHECKING
 IS DONE TO THIS RESPECT.
 */
void computeSHEigMatrixSymmetric( const unsigned int L, SHMatrixType& eig )
{
   // Set the proper size:
   eig.set_size( getNumberOfEvenAssociatedLegendrePolynomials(L), getNumberOfEvenAssociatedLegendrePolynomials(L) );
   // Fill with zeros:
   eig.fill( 0.0f );
   // Set the proper values on the diagonal of the matrix:
   unsigned int pos = 0; // Auxiliar position counter
   for( unsigned int l=0; l<=L/2; ++l ){ // For each degree
      for( int m=-(int)(2*l); m<=(int)(2*l); ++m, ++pos )
         eig(pos,pos) = -(double)( (2*l)*(2*l+1) );
   }
}
   
/** This function generates a diagonal matrix (frt) whose diagonal entries are
 the eigenvalues of the Spherical Harmonics basis functions associated to the
 Funk-Radon transform (FRT) operator (SH are also eigenfunctions with respect
 to this linear operator. The eigenvalues for the SH functions up to degree L
 are computed and stored in the diagonal of eig with the usual order: (l=0,m=0),
 (l=1,m=-1), ..., (l=1,m=1), (l=2,m=-2), ..., (l=2,m=2), ... (l=L,m=-L), 
 ...(l=L,m=L).

 Once again, the corresponding eigenvalue in each case depends only on the
 degree l of the SH basis function and does not depend on the order m. It
 can be computed in terms of the value of Legendre polynomials at x=0:
 
   \lambda_{l,m} = 2·\pi·P_l(0)
 
 NOTE: The eigenvalue associated to odd degrees l is always 0 (the SH
 correspond to antisymmetric functions, and hence the integration in a 
 whole equator yields 0). For even degrees, an alternative expression can
 be obtained:
 
   \lambda_{l,m} = 2·\pi·(-1)^(l/2) (l-1)!!/l!!,
 
 where !! denotes de bouble factorial: (l-1)!! = (l-1)·(l-3)·...·3·1, and
 l!! = l·(l-2)·(l-4)·...·4·2, for even l. Anyway, we made use of the
 recursive computation of Legendre polynomials instead.
 */
void computeSHFRTMatrix( const unsigned int L, SHMatrixType& frt)
{
   // Compute the values of the Legendre polynomials as a block:
   double* buffer = allocateBufferForLegendrePolynomials(L);
   computeLegendrePolynomials( 0.0f, L, buffer );
   // Allocate the memory for the resulting matrix:
   frt.set_size( getNumberOfAssociatedLegendrePolynomials(L), getNumberOfAssociatedLegendrePolynomials(L) );
   // Fill with zeros:
   frt.fill( 0.0f );
   // And place the corresponding values:
   unsigned int pos = 0; // Auxiliar position counter
   for( unsigned int l=0; l<=L; ++l ){ // For each degree
      for( int m=-(int)l; m<=(int)l; ++m, ++pos ){
         frt(pos,pos) = (2*PI) * buffer[l];
      }
   }
   delete[] buffer;
}

/** This function generates a diagonal matrix (frt) whose diagonal entries are
 the eigenvalues of the Spherical Harmonics basis functions associated to the
 Funk-Radon transform (FRT) operator (SH are also eigenfunctions with respect
 to this linear operator. The eigenvalues for the SH functions up to degree L
 are computed and stored in the diagonal of eig with the usual order: (l=0,m=0),
 (l=1,m=-1), ..., (l=1,m=1), (l=2,m=-2), ..., (l=2,m=2), ... (l=L,m=-L),
 ...(l=L,m=L).
 
 Once again, the corresponding eigenvalue in each case depends only on the
 degree l of the SH basis function and does not depend on the order m. It
 can be computed in terms of the value of Legendre polynomials at x=0:
 
   \lambda_{l,m} = 2·\pi·P_l(0)
 NOTE: The eigenvalue associated to odd degrees l is always 0 (the SH
 correspond to antisymmetric functions, and hence the integration in a 
 whole equator yields 0). For even degrees, an alternative expression can
 be obtained:
 
 \lambda_{l,m} = 2·\pi·(-1)^(l/2) (l-1)!!/l!!,
 
 where !! denotes de bouble factorial: (l-1)!! = (l-1)·(l-3)·...·3·1, and
 l!! = l·(l-2)·(l-4)·...·4·2, for even l. Anyway, we made use of the
 recursive computation of Legendre polynomials instead.
 
 As opposed to computeSHFRTMatrix, this function is restricted to even degrees
 of the SH basis, and hence it is only usuful for fucntions defined over the
 unit sphere showing radial symmetry. L IS ASSUME DTO BE EVEN AND NO CHECKING
 IS DONE TO THIS RESPECT.
 */
void computeSHFRTMatrixSymmetric( const unsigned int L, SHMatrixType& frt )
{
   // Compute the values of the Legendre polynomials as a block:
   double* buffer = allocateBufferForEvenLegendrePolynomials(L);
   computeEvenLegendrePolynomials( 0.0f, L, buffer );
   // Allocate the memory for the resulting matrix:
   frt.set_size( getNumberOfEvenAssociatedLegendrePolynomials(L), getNumberOfEvenAssociatedLegendrePolynomials(L) );
   // Fill with zeros:
   frt.fill( 0.0f );
   // And place the corresponding values:
   unsigned int pos = 0; // Auxiliar position counter
   for( unsigned int l=0; l<=L/2; ++l ){ // For each degree
      for( int m=-(int)(2*l); m<=(int)(2*l); ++m, ++pos ){
         frt(pos,pos) = (2*PI) * buffer[l];
      }
   }
   delete[] buffer;
}


/** The exponential integral. Ei(x) is defined as \int_x^\infty \exp(-t)/t dt,
 and hence shows a singularity at x=0. It may be defined in terms of a power
 series expansion (regular part) and the logarithm function (singularity
 at x=0). Although this power series converge as fast as that for the
 exponential function or even faster, it is not an efficient implementation.
 
 Instead, we use the implementation provided by ALGLIB under the GPL license:
 
                  http://www.alglib.net/
 
 This implementation is very similar to that in matlab. Compared to the original
 in ALGLIB, we have stablished much less conservative thresholds, since we
 need a fast performance. Even so, it is accurate enough. If more precission is
 required, the places marked with the word "Tolerance" should be revised.
 
 NOTE: It is assumed that the input, x, is positive. Otherwise, the
 logarithm of x is returned (NaN).
 
 NOTE(2): We have dropped the argument "n" in the original ALGLIB routine, 
 since we need to compute only E_1.
 */
double expint( const double x )
{
   double result;
    double r;
    double t;
    double yk;
    double xk;
    double pk;
    double pkm1;
    double pkm2;
    double qk;
    double qkm1;
    double qkm2;
    double psi;
    double z;
   int k;
    double big;
   unsigned int cont = 0;
    
    big = 1.44115188075855872e17;
    if( x<=1e-3f ){ // asymptotic approximation -- Tolerance
      // Since x is quite small, the Taylor series expansion can be computed with
      // very few coefficients (we keep the first three terms; the modulus
      // of the next term, k=4, will be less than (1e-3)^4/(4·4!) ~ 1e-14:
      return ( -::log(x) - GAMMAEULER + x - x*x/4 + x*x*x/18 );
   }
   if( x>=30.0f ) // For greater values, it is practically zero (below 3e-15) -- Tolerance
      return 0.0f;
   if( x<=1.0f ){
      psi = - GAMMAEULER - ::log(x);
      z   = -x;
        xk  = 0;
        yk  = 1;
        pk  = 0;
        result = 0.0;
        do{
         xk      = xk + 1;
         yk      = yk * z / xk;
         pk     += 1;
         result += yk/pk;
            if( ::fabs(result)>1e-9 ) // Tolerance
            t = fabs(yk/result);
            else
            t = 1;
        }
        while( t>=1e-6 && cont++<100 ); //Tolerance
        result = psi-result;
        return result;
    }
    else{
      k      = 1;
      pkm2   = 1;
        qkm2   = x;
        pkm1   = 1.0;
        qkm1   = x + 1;
        result = pkm1/qkm1;
      do{
         if( ++k%2 == 1 ){
                yk = 1;
                xk = 1 + double(k-1)/double(2);
            }
            else{
            yk = x;
            xk = double(k)/double(2);
         }
         pk = pkm1*yk + pkm2*xk;
         qk = qkm1*yk + qkm2*xk;
            if( ::fabs(qk)>1e-9 ){ // Tolerance
            r = pk/qk;
            t = fabs((result-r)/r);
            result = r;
         }
            else
            t = 1;
         pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            if( fabs(pk)>=big ){
                pkm2 = pkm2/big;
                pkm1 = pkm1/big;
                qkm2 = qkm2/big;
                qkm1 = qkm1/big;
            }
        }
        while( t>=1e-6 && cont++<100 ); //Tolerance
        result = result*exp(-x);
    }
    return result;
}
   
/** The non-singular exponential integral, defined as the power series
 in the definition of the exponential integral without the logarithmic
 singularity and the Euler constant. Although the convergence of this 
 power series expansion is quite fast, it is preferable to use the
 expint routine (for rasonable values) instead.
 
 The relative error in the computation of the Ein function is at most
 1.5e-4 for x=6. The absolut error is nearly 3.6e-4 for x=6.
 
 This errors can be decreased (or incrased, with the advantage
 of faster performance) modifying the thresholds used to decied wether
 to use expint or not.
 */
double Ein( const double x )
{
   if( x<0.1f ){
      // x is very small, so the Taylor series expansion with 3 terms is the most efficient
      // implementation; note that the modulus of the first neglected term (k=4) is at most
      // (1e-1)^4/(4·4!) ~ 1e-6. This means that the maximum error is, at most, in the order
      // of 1e-5.
      return( -x + x*x/4 - x*x*x/18 );
   }
   else if( x>6.0f ){
      // For larger values of x, the expint function is negligible compared to the logarithm,
      // so it makes no sense to compute its value; the relative error with this approximation
      // is at most in the order of 5e-4
      return( -::log(x) - GAMMAEULER );
   }
   else // The normal case. The computation is based on the expint function
      return ( -expint(x) - ::log(x) - GAMMAEULER );
}
   
/** This is quite a specialized function to compute the angular coordinates theta and phi of
 a spherical coordinates system from the cartesian coordinates (x,y,z). Physics convention is
 used, so that the 'z' axis correspond to theta=0, the 'x' axis to (theta=0, phi=0), and the
 'y' axis to (theta=0, phi=PI/2).
 
 Regardless on the actual 'r' coordinate, the vector [x,y,z] is normalized so that it has
 norm 1. This because this function is designed for functions defined over the unit sphere
 (which are those we are able to represent by means of SH expansions).
 
 Besides, this function is specifically designed to work with functions showing radial 
 symmetry (only EVEN degrees of the SH basis functions/associated Legendre polynomials).
 For convenience, all points are translated to the hemisphere defined by 0 <= phi < pi.
 Each point lying in the opposite hemisphere is projected onto its radially symmetric
 point.
 
 NOTE: The tolerances for floating point comparisons are not conservative at all; this
 is because this function is originally designed to translate MRI machinery-provided
 diffusion directions to angular coordinates, and considering a higher precission is 
 n0t necessary at all.
 */
void computeShericalCoordsFromCartesian( const double x, const double y, const double z, double& theta, double& phi )
{
   //===============================================================================================================
   // First of all, eliminate the radial coordinate and work with normalized coordinates:
   double r  = ::sqrt( x*x + y*y + z*z );
   if( r <= 1e-6 ){ // The vector has norm zero
      theta = phi = 0.0f;
      return;
   }
   double x0 = x/r;
   double y0 = y/r;
   double z0 = z/r;
   //===============================================================================================================
   // The computation of the theta coordinate is quite simple, since it only depends on z0:
   theta     = ::acos( z0 );
   // If the theta coordinate corresponds to the 'z' or '-z' axis, the phi coordinate is undefined:
   if(   ( theta<1e-3 )   ||   ( PI-theta<1e-3 )   )
      phi = 0.0f;
   else{
      //===============================================================================================================
      // If the direction is not parallel to the 'z' axis, the phi coordinate has to be actually computed. Since, with
      // phisics convention and for radius 1, x0 = sin(theta)·cos(phi), y0 = sin(theta)·sin(phi), we can divide by
      // sin(theta) to obtain the cosine and the sine of phi. Since theta is not 0 or pi, the division is well defined.
      x0       /= ::sin( theta );
      y0       /= ::sin( theta );
      //===============================================================================================================
      // To avoid numerical issues, we first check if the corresponding phi coordinante correspond to a direction
      // aligned with either the 'x' or the 'y' axis.
      if( ::fabs(x0)<0.01 ){
         if( y0>0.0f )
            phi = PI/2;          // 'y'  axis
         else
            phi = 3*PI/2;        // '-y' axis
      }
      else if( ::fabs(y0)<0.01 ){
         if( x0>0.0f )
            phi = 0.0f;          // 'x'  axis
         else
            phi = PI;            // '-x' axis
      }
      else{ // No numerical issues. Proceed normally
         //===============================================================================================================
         // The value of phi can be obtained from x0 with the acos or from y0 with the asin. To increase the accuacy in
         // the computation, both values are computed and averaged:
         double p1 = ::acos( x0 ); //   0   <= p1 <=  PI
         double p2 = ::asin( y0 ); // -PI/2 <= p2 <= PI/2
         // However, we have to avoid the ambiguity in the computation of the inverse sine and cosine (each value between
         // -1 and 1 may correspond to the sine/cosine of two different angles in the range [0,2*pi]):
         if( x0<0.0f ){
            if( y0<0.0f ){ // Third quadrant
               p1 = 2*PI - p1; // The angle in the second quadrant is returned by acos; reflex the angle with respect to p1=PI
               p2 = PI   - p2; // The angle in the fourth quadrant is returned by asin; reflex with respect to p2=3PI/2
            }
            else           // Second quadrant
               p2 = PI   - p2; // The angle in the first quadrant is returned by asin. Compute the supplementary angle
         }
         else{
            if( y0<0.0f ){ // Fourth quadrant
               p1 = 2*PI - p1; // The angle in the first quadrant is returned by acos; invert p1 and place in the range [0,2PI]
               p2 = 2*PI + p2; // The angle in the fourth quadrant is returned; place in the range [0,2PI]
            }
         }
         // Average p1 and p2 to compute phi:
         phi = (p1+p2)/2;
      } // else [ if( ::fabs(x0)<0.02 ) ]
   } // else [  if(   ( theta<1e-3 )   ||   ( PI-theta<1e-3 )   )  ]
   //===============================================================================================================
   // The last step is to project the points int the hemisphere phi>PI to points in the hemispher phi<PI. The '-z'
   // axis is also projected onto its radially symmetric position, 'z'.
   if( (PI-theta)<1e-3 ) // This is the '-z' axis
      theta = phi = 0.0f;
   else if( phi<0.0f )   // This situation should never be reached...
      phi   = 0.0f;
   else if( phi>=PI ){   // Hemisphere with phi>pi
      theta = PI - theta;
      phi   = phi - PI;
   }
}
   
/** This function sequentially uses:
 
 ComputeShericalCoordsFromCartesian( const double x, const double y, const double z, double& theta, double& phi )
 
 to obtain the angular sphercial coordinates theta and phi of a set of points given in Cartesian coordinates
 by the vectors x, y, and z with length N (physics conventions). See the documentation above for details.
 */
 
void computeShericalCoordsFromCartesian( const double* x, const double* y, const double* z, double* theta, double* phi, const unsigned int N )
{
   for( unsigned int n=0; n<N; ++n )
      computeShericalCoordsFromCartesian( x[n], y[n], z[n], theta[n], phi[n] );
}

   
} // End namespace shmaths

#endif // #ifndef _sphericalHarmonics_cxx
