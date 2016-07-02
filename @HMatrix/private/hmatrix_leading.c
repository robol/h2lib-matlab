#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));

  /* Get the leading block of A, and the associated cluster trees. */
  phmatrix Al = NULL;
  ref_hmatrix (&Al, A->son[0]);
  
  pccluster lrc = A->rc->son[0];
  pccluster lcc = A->cc->son[0];

  /* Copy the relevant data in the output */
  SERIALIZE_POINTER_TO_PROPERTY(prhs[1], "hmatrix", Al);
  SERIALIZE_POINTER_TO_PROPERTY(prhs[2], "cluster", lrc);
  SERIALIZE_POINTER_TO_PROPERTY(prhs[3], "cluster", lcc);

} 
