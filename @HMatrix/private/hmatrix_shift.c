#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  phmatrix Ar = clone_hmatrix (A);
  
  field r = *mxGetPr(prhs[1]);
  
#ifdef USE_COMPLEX
  double * ri = mxGetPi(prhs[1]);
  if (ri)
    r = r + I * *ri;
#endif
  
  shift_hmatrix (r, Ar);
  
  SERIALIZE_POINTER (plhs[0], Ar);
}
