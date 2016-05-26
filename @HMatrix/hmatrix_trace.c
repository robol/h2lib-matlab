#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray* plhs[], 
		 int nrhs, const mxArray* prhs[])
{
  pchmatrix H = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  plhs[0] = mxCreateDoubleScalar (hmatrix_trace(H));
}
