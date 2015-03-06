#include "mex.h"
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <assert.h>
#include "hmatlab.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int n = mxGetM(prhs[0]);
  int* a = malloc(sizeof(int)*n);
  int i;

  for(i = 0;i < n;i++)
    a[i]=i;

  pcluster rc = create_cluster(a,n);
  pcluster cc = create_cluster(a,n);

  update_cluster(rc);
  update_cluster(cc);

  phmatrix A = create_tridiag_hmatrix (mxGetPr(prhs[0]), mxGetPr(prhs[1]), 
			               mxGetPr(prhs[2]), rc, cc);

  if (nlhs > 0)
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	*((long long *) mxGetData(plhs[0])) = (long long) A;
    }
}
