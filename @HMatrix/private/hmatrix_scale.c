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
  phmatrix Ascaled = clone_hmatrix (A);
  
  double * sr = mxGetPr(prhs[1]);
  double * si = mxGetPi(prhs[1]);

#ifdef USE_COMPLEX
  scale_hmatrix (*sr + I * (si ? *si : 0), Ascaled);
#else
  scale_hmatrix (*sr, Ascaled);
#endif
  

  SERIALIZE_POINTER (plhs[0], Ascaled);
}
