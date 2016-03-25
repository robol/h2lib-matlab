#include "mex.h"
#include <hmatrix.h>
#include <avector.h>
#include <string.h>
#include <hmatlab.h>

#ifndef false
#define false 0
#endif

/*
 * Compute the matrix product A v = y. 
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int i;
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  
  /* Real and imaginary part of the input vector v */
  double * vr = mxGetPr(prhs[1]);
  double * vi = mxGetPi(prhs[1]);
  
  int n = mxGetM(prhs[1]);
  
  /* Create the space to hold the input data in the H2Lib 
   * calls. We need both input and output. */
  field * v = mxMalloc (sizeof (field) * n);
  field * y = mxMalloc (sizeof (field) * n);
  
  /* We take a copy of the input vector */
  for (i = 0; i < n; i++) {
#ifdef USE_COMPLEX  
    v[i] = vr[i] + I * (vi ? vi[i] : 0);
#else
	v[i] = vr[i];
#endif
  }

  avector av;
  pavector pav = init_pointer_avector(&av, v, n);

  avector ay;  
  pavector pay = init_pointer_avector (&ay, y, n);
  memset (y, 0, sizeof (field) * n);

  mvm_hmatrix_avector (1.0, false, A, pav, pay);

#ifdef USE_COMPLEX
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);  
#else
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);  
#endif
  
  /* Copy the result back in the output vector */
  double * yr = mxGetPr(plhs[0]);
  double * yi = mxGetPi(plhs[0]);
  
  if (yi)
      memset (yi, 0, sizeof (double) * n);
  
  for (i = 0; i < n; i++) {
#ifdef USE_COMPLEX
    yr[i] = creal(y[i]);
    yi[i] = cimag(y[i]);  
#else
    yr[i] = y[i];
#endif 
  }
  
  mxFree (v);
  mxFree (y);

  uninit_avector(pav);
  uninit_avector(pay);
}
