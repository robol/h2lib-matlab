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
  phmatrix B = DESERIALIZE_POINTER (mxGetProperty (prhs[1], 0, "hmatrix"));

  ptruncmode tm = new_releucl_truncmode();

  phmatrix C = NULL; ref_hmatrix(&C,clone_hmatrix (A));

  add_hmatrix(- FIELD_ONE, B, tm, h2lib_eps, C);

  SERIALIZE_POINTER (plhs[0], C);
}
