#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>


#define eps 0.00000001

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = (phmatrix) *((phmatrix*) mxGetData(prhs[0]));
  phmatrix B = (phmatrix) *((phmatrix*) mxGetData(prhs[1]));
  ptruncmode  	tm;
  tm=new_releucl_truncmode();
  add_hmatrix(1.0,A,tm,eps,B);
  plhs[0] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
  *((long long *) mxGetData(plhs[0])) = (long long) B;
}
