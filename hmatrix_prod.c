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
  phmatrix B = (phmatrix) *((phmatrix*) mxGetData(prhs[1]));
  phmatrix C = NULL;
  double * a = malloc (A->rc->size * sizeof(double));
  memset (a, 0, sizeof(double) * A->rc->size);
  C = create_tridiag_hmatrix (a, a, a, A->rc, A->cc);
  free (a);
  /* C = new_hmatrix(A->rc,A->cc); */
  ptruncmode  tm;
  tm=new_releucl_truncmode();
  addmul_hmatrix(1.0,false,A,false,B,tm,eps,C);
  plhs[0] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
  *((long long *) mxGetData(plhs[0])) = (long long) C;
}
