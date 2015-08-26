#include "mex.h"
#include <h2matrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  ph2matrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "h2matrix"));

  if (nrhs == 1)
    {
      plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
      double * size = mxGetData(plhs[0]);

      size[0] = A->rb->t->size;
      size[1] = A->cb->t->size;
    }
  else
    {
      int dimension = *mxGetPr(prhs[1]);

      switch (dimension) 
	{
	case 1:
	  plhs[0] = mxCreateDoubleScalar (A->rb->t->size);
	case 2:
	  plhs[0] = mxCreateDoubleScalar (A->cb->t->size);
	default:
	  plhs[0] = mxCreateDoubleScalar (1);
	}
    }
}
