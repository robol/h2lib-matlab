#include "mex.h"
#include <hmatrix.h>
#include <avector.h>
#include <string.h>
#include <hmatlab.h>

#ifndef false
#define false 0
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  /* phmatrix A = (phmatrix) *((phmatrix*) mxGetData(prhs[0])); */
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  double * v = mxGetPr(prhs[1]);

  int n = mxGetM(prhs[1]);

  avector av;
  pavector pav = init_avector(&av, n);
  pav->v = v;

  avector y;
  pavector py = init_avector (&y, n);
  fill_avector (py, 0);

  mvm_hmatrix_avector (1.0, false, A, pav, py);

  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  memcpy (mxGetPr(plhs[0]), py->v, n * sizeof(double));  
}
