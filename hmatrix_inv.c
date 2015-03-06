#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include "hmatlab.h"

#define eps 0.00000001

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = (phmatrix) *((phmatrix*) mxGetData(prhs[0]));
  double * a = malloc (sizeof(double) * A->rc->size);
  memset (a, 0, sizeof(double) * A->rc->size);
  phmatrix B = create_tridiag_hmatrix(a,a,a,A->rc,A->cc);
  free (a);
  ptruncmode  tm;
  tm=new_releucl_truncmode();
  invert_hmatrix(A,B,tm,eps);

  del_hmatrix (B);
  plhs[0] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
  *((long long *) mxGetData(plhs[0])) = (long long) A;
}
