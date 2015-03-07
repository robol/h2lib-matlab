#include "mex.h"
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  /* Note: Here  prhs[0] is the cluster object itself. */
  int n = *mxGetPr(prhs[1]);
  int k = *mxGetPr(prhs[2]);

  int * a = malloc (sizeof (int) * n);
  int i;

  for (i = 0; i < n; i++)
    a[i] = i;

  pcluster cluster = std_subdivision_scheme_cluster (a, n, k);

  SERIALIZE_POINTER_TO_PROPERTY ((mxArray*) prhs[0], "cluster", cluster);
}
