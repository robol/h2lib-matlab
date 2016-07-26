#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include <hmatlab.h>
#include <harith.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     int n = mxGetN(prhs[1]);
     int m = mxGetM(prhs[1]);
     int i;
     int dims[] = { m, n };

     mexPrintf("m = %d, n = %d\n", m, n);

     phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
     amatrix am, sm;
     
     ptruncmode tm = new_releucl_truncmode();
     phmatrix Alr=clone_hmatrix(A);

     /* Prepare the output with the same shape of the input */
#ifdef USE_COMPLEX
     plhs[0] = mxCreateDoubleMatrix (m, n, mxCOMPLEX);
#else
     plhs[0] = mxCreateDoubleMatrix (m, n, mxREAL);
#endif

     double * sr = mxGetPr(prhs[1]);
     double * si = mxGetPi(prhs[1]);

     double * or = mxGetPr(plhs[0]);
     double * oi = mxGetPi(plhs[0]);

#ifndef USE_COMPLEX
     init_pointer_amatrix(&sm, sr, m, n);
     init_pointer_amatrix(&am, or, m, n);
#else
     init_amatrix(&sm, m, n);        
     init_amatrix(&am, m, n);
     for (i = 0; i < m * n; i++) {
	 sm.a[i] = sr[i] + I * (si ? si[i] : 0.0);     
     }
#endif

     memset (am.a, 0, sizeof(field) * m * n);
     addmul_hmatrix_amatrix_amatrix(1.0, false, A, false, &sm, false, &am);

#ifdef USE_COMPLEX
     for (i = 0; i < m * n; i++) {
	 or[i] = creal(am.a[i]);
	 oi[i] = cimag(am.a[i]);
     }
#endif

     uninit_amatrix(&am);
     uninit_amatrix(&sm);

}
