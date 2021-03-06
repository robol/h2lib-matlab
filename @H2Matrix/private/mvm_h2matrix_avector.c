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
  ph2matrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "h2matrix"));
  double * v = mxGetPr(prhs[1]);

  int n = mxGetM(prhs[1]);

  avector av;
  pavector pav = init_pointer_avector(&av, v, n);

  avector y;
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  pavector py = init_pointer_avector (&y, mxGetPr(plhs[0]), n);
  fill_avector (py, 0);

  mvm_h2matrix_avector (1.0, false, A, pav, py);

  uninit_avector(pav);
  uninit_avector(py);
}
