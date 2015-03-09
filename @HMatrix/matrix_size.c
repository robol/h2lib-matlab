#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));

  if (nrhs == 1)
    {
      plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
      double * size = mxGetData(plhs[0]);

      size[0] = A->rc->size;
      size[1] = A->cc->size;
    }
  else
    {
      int dimension = *mxGetPr(prhs[1]);

      switch (dimension) 
	{
	case 1:
	  plhs[0] = mxCreateDoubleScalar (A->rc->size);
	case 2:
	  plhs[0] = mxCreateDoubleScalar (A->cc->size);
	default:
	  plhs[0] = mxCreateDoubleScalar (1);
	}
    }
}
