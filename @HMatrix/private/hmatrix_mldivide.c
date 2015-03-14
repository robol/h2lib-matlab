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
     int i;
     int dims[] = { 1, n };

     phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));

     ptruncmode tm = new_releucl_truncmode();
     phmatrix Alr=clone_hmatrix(A);
     lrdecomp_hmatrix(Alr,tm,h2lib_eps);

     plhs[0] = mxCreateNumericArray (2, dims, mxINT64_CLASS, mxREAL);
     long * pointers = (long*) mxGetData(plhs[0]);

     for (i = 0; i < n; i++)
       {
         phmatrix B = DESERIALIZE_POINTER (mxGetProperty (prhs[1], i, "hmatrix"));
         phmatrix Blr=clone_hmatrix(B);
         triangularinvmul_hmatrix(true,true,false,Alr,tm,h2lib_eps,false,Blr);
         triangularinvmul_hmatrix(false,false,false,Alr,tm,h2lib_eps,false,Blr);
         pointers[i] = (long) Blr;
       }

     del_hmatrix(Alr);
}
