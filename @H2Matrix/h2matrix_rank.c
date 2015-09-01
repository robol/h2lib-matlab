#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray* plhs[], 
		 int nrhs, const mxArray* prhs[])
{
  pch2matrix H = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "h2matrix"));
  plhs[0] = mxCreateDoubleScalar (h2matrix_get_rank (H));
}
