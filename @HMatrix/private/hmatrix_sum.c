#define USE_COMPLEX 1

#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  phmatrix B = DESERIALIZE_POINTER (mxGetProperty (prhs[1], 0, "hmatrix"));
  ptruncmode tm = new_releucl_truncmode();

  phmatrix C = NULL;  ref_hmatrix(&C, clone_hmatrix( B) );
 
  real eps = h2lib_eps;

  /* printf ("Passing eps = %e\n", h2lib_eps); */
  add_hmatrix(FIELD_ONE, A, tm, eps, C);

  SERIALIZE_POINTER (plhs[0], C);
}
