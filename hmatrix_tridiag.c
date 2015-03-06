#include "mex.h"
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <assert.h>

#define MIN_SIZE_THRESHOLD 4
#define MATRIX_ELEM(A,i,j,n) A[(i)*(n)+(j)]

static pcluster create_cluster(int *a, int n)
{
  pcluster clust=NULL;

  if(n <= MIN_SIZE_THRESHOLD)
    {
      clust = new_cluster(n, a, 0, 1);
    }
  else
    {
      cluster * c = NULL;
      clust=new_cluster(n,a,2,1);
      c = clust;
      c->son[0]=create_cluster(a,n/2);
      c->son[1]=create_cluster(a+n/2,n-n/2);
  }

  return clust;
}

phmatrix create_hmatrix (double * a, double * b, double * c, 
			 pcluster rc, pcluster cc)
{
  int n = rc->size;
  phmatrix A = NULL;

  assert (rc->size == cc->size);

  if (n <= MIN_SIZE_THRESHOLD)
    {
      A = new_full_hmatrix (rc, cc);
      int i;

      memset (A->f->a, 0, sizeof (double) * n * n);

      A->f->a[0] = a[0];
      for (i = 0; i < n - 1; i++)
	{
	  MATRIX_ELEM(A->f->a, i+1, i, n) = b[i];
	  MATRIX_ELEM(A->f->a, i + 1, i + 1, n) = a[i+1];
	  MATRIX_ELEM(A->f->a, i, i + 1, n) = c[i];
	}

      update_hmatrix (A);
    }
  else 
    {
      /* Perform subdivision */
      A = new_super_hmatrix (rc, cc, 2, 2);

      phmatrix A11 = create_hmatrix (a, b, c, rc->son[0], cc->son[0]);
      phmatrix A22 = create_hmatrix (a + rc->son[0]->size, 
				     b + rc->son[0]->size,
				     c + rc->son[0]->size,
				     rc->son[1], cc->son[1]);

      phmatrix A12 = new_rk_hmatrix (rc->son[0], cc->son[1], 1);
      phmatrix A21 = new_rk_hmatrix (rc->son[1], cc->son[0], 1);

      memset (A12->r->A.a, 0, sizeof (double) * rc->son[0]->size);
      memset (A12->r->B.a, 0, sizeof (double) * cc->son[1]->size);
      memset (A21->r->A.a, 0, sizeof (double) * rc->son[1]->size);
      memset (A21->r->B.a, 0, sizeof (double) * cc->son[0]->size);

      /* Fill in the low rank parts */
      A12->r->A.a[rc->son[0]->size - 1] = b[cc->son[0]->size - 1];
      A12->r->B.a[0] = 1.0;
      A21->r->A.a[0] = c[rc->son[0]->size - 1];
      A21->r->B.a[cc->son[0]->size - 1] = 1.0;

      A->son[0] = A11;
      A->son[1] = A21;
      A->son[2] = A12;
      A->son[3] = A22;

      update_hmatrix (A);
    }

  return A;
}

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

  phmatrix A = create_hmatrix (mxGetPr(prhs[0]), mxGetPr(prhs[1]), 
			       mxGetPr(prhs[2]), rc, cc);

  if (nlhs > 0)
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	*((long long *) mxGetData(plhs[0])) = (long long) A;
    }
}
