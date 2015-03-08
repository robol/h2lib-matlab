#include "mex.h"
#include <cluster.h>
#include <hmatlab.h>

void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrsh, const mxArray* prhs[])
{
  pcluster c = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "cluster"));

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double * size = mxGetPr(plhs[0]);
  *size = c->size;
}
