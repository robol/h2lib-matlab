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

     /* Prepare the output with the same shape of the input */
     plhs[0] = mxCreateDoubleMatrix (m, n, mxREAL);

     /* Copy the input into the output, since it will be overwritten. */
     memcpy (mxGetPr(plhs[0]), mxGetPr(prhs[1]), n * m * sizeof (double));

     am.rows = m;
     am.cols = n;
     am.ld   = m;
     am.owner = NULL;
     am.a = mxGetPr(plhs[0]);

     triangularinvmul_hmatrix_amatrix (true, true, false, Alr, false, &am);
     triangularinvmul_hmatrix_amatrix (false, false, false, Alr, false, &am);

     del_hmatrix(Alr);
}
