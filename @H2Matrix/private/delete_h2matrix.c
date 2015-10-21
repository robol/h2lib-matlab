#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[])
{
  ph2matrix H = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "h2matrix"));
  del_h2matrix (H);
}
