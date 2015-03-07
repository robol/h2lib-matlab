#include "mex.h"
#include <hmatrix.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  plhs[0] = mxCreateNumericMatrix(1, 2, mxINT64_CLASS, mxREAL);
  long * size = mxGetData(plhs[0]);
  phmatrix A = (phmatrix) *((phmatrix*) mxGetData(prhs[0]));

  size[0] = A->rc->size;
  size[1] = A->cc->size;
}
