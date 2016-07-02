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
     phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));

     phmatrix L = NULL; 
     ref_hmatrix(&L, clone_hmatrix(A));
     phmatrix U = NULL;
     ref_hmatrix(&U, clone_hmatrix(A));

     identity_hmatrix(L);
     identity_hmatrix(U);
     ptruncmode tm = new_releucl_truncmode();
     phmatrix Alr=clone_hmatrix(A);
     lrdecomp_hmatrix(Alr,tm,h2lib_eps);
     triangularmul_hmatrix 	(true,true,false,Alr,tm,h2lib_eps,false,L);
     triangularmul_hmatrix 	(false,false,false,Alr,tm,h2lib_eps,false,U);
     del_hmatrix(Alr);

     SERIALIZE_POINTER (plhs[0], L);
     SERIALIZE_POINTER (plhs[1], U);
}
