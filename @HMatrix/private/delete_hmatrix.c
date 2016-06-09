#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[])
{
  phmatrix H = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));

  unref_hmatrix(H);
}
