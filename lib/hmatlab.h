#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>

#ifndef _HMATLAB_H
#define _HMATLAB_H

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
 * @param b The superdiagonal entries of the matrix.
 * @param c The subdiagonal entries of the matrix. 
 */
phmatrix create_tridiag_hmatrix (double * a, double * b, double * c, 
				 pccluster rc, pccluster cc);

#define SERIALIZE_POINTER_TO_PROPERTY(a, property, ptr) {		\
    mxArray * p = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);	\
    *((long *) mxGetData(p)) = ((long) (ptr));				\
    mxSetProperty((mxArray*) a, 0, property, p);				\
  }


#define SERIALIZE_POINTER(a, ptr) {					\
    mxArray ** b = &(a);						\
    *b = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);		\
    *((long *) mxGetData(*b)) = ((long) (ptr));				\
  }

#define DESERIALIZE_POINTER(a) (void*) *((long *) mxGetData(a))

#endif
