#include "mex.h"
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <assert.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  pcluster rc = DESERIALIZE_POINTER (mxGetProperty (prhs[1], 0, "cluster"));
  pcluster cc = DESERIALIZE_POINTER (mxGetProperty (prhs[2], 0, "cluster"));

  phmatrix A = create_band_hmatrix (mxGetPr(prhs[3]), mxGetPr(prhs[4]), 
			               mxGetPr(prhs[5]), *mxGetPr(prhs[6]), *mxGetPr(prhs[7]), rc, cc);

  SERIALIZE_POINTER_TO_PROPERTY (prhs[0], "hmatrix", A);
}
