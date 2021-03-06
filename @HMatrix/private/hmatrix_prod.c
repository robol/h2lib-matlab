#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include "hmatlab.h"

#ifndef false
#define false 0
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty(prhs[0], 0, "hmatrix"));
  phmatrix B = DESERIALIZE_POINTER (mxGetProperty(prhs[1], 0, "hmatrix"));

  ptruncmode  tm = new_releucl_truncmode();
  phmatrix C = NULL;

  double * a = malloc (A->rc->size * sizeof(double));
  memset (a, 0, sizeof(double) * A->rc->size);

  ref_hmatrix(&C,  create_real_tridiag_hmatrix (a, a, a, A->rc, A->cc)); 
  free (a);

  addmul_hmatrix(FIELD_ONE, false, A, false, B, tm, h2lib_eps, C);

  SERIALIZE_POINTER (plhs[0], C);
}
