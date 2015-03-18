#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>

#ifndef _HMATLAB_H
#define _HMATLAB_H

#define h2lib_eps 1e-13

/**
 * @brief Create an index cluster representing the standard subdivision
 * scheme obtained by recursively subdividing the segment by bisection
 * until a minimum length is reached.
 *
 * @param a The array containing the indices.
 * @param n The length of the array <code>a</code>.
 * @param k The maximum admissibile length for a segment to not
 * perform further bisection.
 *
 * @return A newly allocated cluster.
 */
pcluster std_subdivision_scheme_cluster (int * a, int n, int k);

/**
 * @brief Create a hierarichal matrix representing a tridiagonal
 * matrix.
 *
 * @param a The diagonal entries of the matrix.
 * @param b The subdiagonal entries of the matrix.
 * @param c The superdiagonal entries of the matrix.
 */
phmatrix create_tridiag_hmatrix (double * a, double * b, double * c,
				 pccluster rc, pccluster cc);

/**
 * @brief Create a hierarichal matrix representing a (p,q)-band matrix
 * 
 *
 * @param a The diagonal entries of the matrix.
 * @param b The subdiagonal entries of the matrix. Stored as the subcolumns of the matrix and with the last n-p columns filled by zeros to allow recursion.
 * @param c The superdiagonal entries of the matrix. Stored as the superrows of the matrix and with the last n-q rows filled by zeros to allow recursion.
 */

phmatrix create_band_hmatrix (double * a, double * b, double * c, int p, int q,
			         pccluster rc, pccluster cc);

/**
 * @brief Create a hierarichal matrix representing a (U^t,V^t,W^t,Z^t)-generator matrix ( i.e. tril(A)= U^t V, triu(A)=W^t Z and diag(A)=d )
 * 
 *
 * @param d The diagonal entries of the matrix.
 * @param U,V the subdiagonal generators  of the matrix. 
 * @param W,Z the superdiagonal generators  of the matrix.
 * @param ksub,ksup the sub and superdiagonal rank of the matrix.
 */

phmatrix create_generators_hmatrix (double * d, double * U, double * V, double * W, double * Z,
			         int ksub, int ksup, pccluster rc, pccluster cc);

/**
 * @brief Obtain the quasiseparable rank of the hierarchical matrix H. 
 * 
 * Currently this routine is implemented in a suboptimal way so that it
 * returns the rank of the first son of H. This can be a good guess for
 * quasiseparability rank in the general case, but there is guarantee. 
 *
 * A good implementation should walk all the tree and take the maximum on the
 * ranks. 
 * 
 * @param H The hmatrix whose rank should be obtained. 
 * @return The QS rank of H.
 */
size_t hmatrix_get_rank (pchmatrix H);

/**
 * @brief Obtain the number of rows of H. 
 *
 * @param H A pointer to the hmatrix struct.
 * @return An integer representing the number of rows. 
 */
size_t hmatrix_get_m (pchmatrix H);

/**
 * @brief Obtain the number of columns of H. 
 *
 * @param H A pointer to the hmatrix struct.
 * @return An integer representing the number of columns. 
 */
size_t hmatrix_get_n (pchmatrix H);

/**
 * @brief Macro used to serialize pointers in classdef-objects properties
 * in MATLAB. 
 */
#define SERIALIZE_POINTER_TO_PROPERTY(a, property, ptr) {		\
    mxArray * p = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);	\
    *((long *) mxGetData(p)) = ((long) (ptr));				\
    mxSetProperty((mxArray*) a, 0, property, p);				\
  }

/**
 * @brief Macro used to serialize pointers to mxArrays in MATLAB. 
 */
#define SERIALIZE_POINTER(a, ptr) {					\
    mxArray ** b = &(a);						\
    *b = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);		\
    *((long *) mxGetData(*b)) = ((long) (ptr));				\
  }

/**
 * @brief Marco used to deserialize pointers from mxArray in MATLAB. 
 */
#define DESERIALIZE_POINTER(a) (void*) *((long *) mxGetData(a))

#endif
