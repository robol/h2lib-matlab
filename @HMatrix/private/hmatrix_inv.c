#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include "hmatlab.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  phmatrix Ainv = NULL;
  ref_hmatrix(&Ainv ,clone_hmatrix (A));

  /* Allocate space for the inversion */
  double * a = malloc (sizeof(double) * A->rc->size);
  memset (a, 0, sizeof(double) * A->rc->size);
  phmatrix B = NULL;
  ref_hmatrix(&B, create_real_tridiag_hmatrix(a,a,a,A->rc,A->cc));
  free (a);

  ptruncmode  tm;
  tm=new_releucl_truncmode();

  invert_hmatrix (Ainv, B, tm, h2lib_eps);

  unref_hmatrix (B);

  SERIALIZE_POINTER (plhs[0], Ainv);
}
