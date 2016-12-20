#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix H = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  mwSize m, n, k, i;

  prkmatrix NE = H->son[1]->r;
  k = NE->k;
  m = NE->A.rows;
  n = NE->B.cols;

#ifdef USE_COMPLEX
  plhs[0] = mxCreateDoubleMatrix(m, k, mxCOMPLEX);
  plhs[1] = mxCreateDoubleMatrix(n, k, mxCOMPLEX);

  double *Ar = mxGetPr(plhs[0]), *Ai = mxGetPi(plhs[0]);
  double *Br = mxGetPr(plhs[1]), *Bi = mxGetPi(plhs[1]);

  double *A = (double*) NE->A.a;
  double *B = (double*) NE->B.a;

  for (i = 0; i < m * k; i++) {
    Ar[i] = A[2*i];
    Ai[i] = A[2*i+1];
  }

  for (i = 0; i < n * k; i++) {
    Br[i] = B[2*i];
    Bi[i] = B[2*i+1];
  }
    
#else
  plhs[0] = mxCreateDoubleMatrix(m, k, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(n, k, mxREAL);

  memcpy (mxGetPr(plhs[0]), NE->A.a, m * k * sizeof (field));
  memcpy (mxGetPr(plhs[1]), NE->B.a, n * k * sizeof (field));
#endif

} 
