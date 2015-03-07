#include "mex.h"
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int n = 8, i;

  uint * indices = malloc (sizeof (uint) * n);

  for (i = 0; i < n; i++)
    indices[i] = i;

  pcluster rowCluster = new_cluster (n, indices, 2, 1);
  pcluster colCluster = new_cluster (n, indices, 2, 1);

  for (i = 0; i < 2; i++)
    {
      cluster * rc = rowCluster;
      rc->son[0] = new_cluster (n / 2, indices, 2, 1);
      rc->son[1] = new_cluster (n / 2, indices + n/2, 2, 1);

      cluster * cc = colCluster;
      cc->son[0] = new_cluster (n / 2, indices, 2, 1);
      cc->son[1] = new_cluster (n / 2, indices + n/2, 2, 1);
    }

  update_cluster (rowCluster);
  update_cluster (colCluster);

  phmatrix A = new_super_hmatrix(rowCluster, colCluster, 2, 2);

  phmatrix A11 = new_full_hmatrix(rowCluster->son[0], colCluster->son[0]);
  phmatrix A22 = new_full_hmatrix(rowCluster->son[1], colCluster->son[1]);
  phmatrix A12 = new_rk_hmatrix(rowCluster->son[0], colCluster->son[1], 1);
  phmatrix A21 = new_rk_hmatrix(rowCluster->son[1], colCluster->son[0], 1);

  A->son[0] = A11;
  A->son[1] = A21;
  A->son[2] = A12;
  A->son[3] = A22;

  /* Fill the matrices */
  init_amatrix (A11->f, n/2, n/2);
  memcpy (A11->f->a, mxGetPr(prhs[0]), n*n / 4 * sizeof(double));
  init_amatrix (A22->f, n/2, n/2);
  memcpy (A22->f->a, mxGetPr(prhs[1]), n*n / 4 * sizeof(double));

  update_hmatrix (A);

  plhs[0] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
  *((long long *) mxGetData(plhs[0])) = (long long) A;
}
