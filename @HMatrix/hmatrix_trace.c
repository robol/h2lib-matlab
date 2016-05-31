#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray* plhs[], 
		 int nrhs, const mxArray* prhs[])
{
    pchmatrix H = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
#ifdef USE_COMPLEX
  field t = hmatrix_trace(H);
  plhs[0] = mxCreateDoubleMatrix (1, 1, mxCOMPLEX);
  *mxGetPr(plhs[0]) = creal(t);
  *mxGetPi(plhs[0]) = cimag(t);
#else
  plhs[0] = mxCreateDoubleScalar (hmatrix_trace(H));
#endif

}
