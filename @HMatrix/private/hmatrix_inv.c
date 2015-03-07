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
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  phmatrix Ainv = clonestructure_hmatrix (A);
  copy_hmatrix (A, Ainv);

  /* Allocate space for the inversion */
  double * a = malloc (sizeof(double) * A->rc->size);
  memset (a, 0, sizeof(double) * A->rc->size);
  phmatrix B = create_tridiag_hmatrix(a,a,a,A->rc,A->cc);
  free (a);

  ptruncmode  tm;
  tm=new_releucl_truncmode();

  invert_hmatrix (Ainv, B, tm, eps);

  del_hmatrix (B);

  SERIALIZE_POINTER (plhs[0], Ainv);
}