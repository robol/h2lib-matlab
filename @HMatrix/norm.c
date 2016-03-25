#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  
  if (nrhs == 1 || nrhs == 2 && *mxGetPr(prhs[1]) == 2)
      plhs[0] = mxCreateDoubleScalar (norm2_hmatrix (A));
  else
      mexErrMsgTxt("Only the 2-norm is supported");
}
