#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include <hmatlab.h>
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     int n = mxGetN(prhs[1]);
     int m = mxGetM(prhs[1]);
     int i;
     int dims[] = { m, n };

     phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
     amatrix am;
     ptruncmode tm = new_releucl_truncmode();
     phmatrix Alr=clone_hmatrix(A);

     lrdecomp_hmatrix(Alr,tm,h2lib_eps);
     
     double * mr = mxGetPr(prhs[1]);
     double * mi = mxGetPi(prhs[1]);

     /* Prepare the output with the same shape of the input */
#ifdef USE_COMPLEX
     plhs[0] = mxCreateDoubleMatrix (m, n, mxCOMPLEX);     
#else
     plhs[0] = mxCreateDoubleMatrix (m, n, mxREAL);
#endif

     double * or = mxGetPr(plhs[0]);
     double * oi = mxGetPi(plhs[0]);

     am.rows = m;
     am.cols = n;
     am.ld   = m;
     am.owner = NULL;     
     
     /* Take a copy of the result, so we have it in the correct HMatrix format. */
	am.a = mxMalloc (sizeof (field) * n * m);
	
	for (i = 0; i < m * n; i++)
#ifdef USE_COMPLEX	
	    am.a[i] = mr[i] + I * (mi ? mi[i] : 0.0);     
#else
  	    am.a[i] = mr[i];     
#endif

     triangularinvmul_hmatrix_amatrix (true, true, false, Alr, false, &am);
     triangularinvmul_hmatrix_amatrix (false, false, false, Alr, false, &am);
     
     /* Copy the result back in place. */
     for (i = 0; i < n * m; i++) {
#ifdef USE_COMPLEX
		or[i] = creal (am.a[i]);
		oi[i] = cimag (am.a[i]);
#else
	    or[i] = am.a[i];		
#endif     
	 }
	 
	 mxFree (am.a);
     del_hmatrix(Alr);
}
