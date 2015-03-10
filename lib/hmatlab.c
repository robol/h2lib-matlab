#include "hmatlab.h"

#define MATRIX_ELEM(A,i,j,n) A[(i)*(n)+(j)]

pcluster std_subdivision_scheme_cluster(int *a, int n, int k)
{
  pcluster clust=NULL;

  if (n <= k)
    {
      clust = new_cluster(n, a, 0, 1);
    }
  else
    {
      clust=new_cluster(n,a,2,1);
      clust->son[0] = std_subdivision_scheme_cluster (a, n/2, k);
      clust->son[1] = std_subdivision_scheme_cluster (a+n/2,n-n/2, k);
  }

  update_cluster (clust);
  return clust;
}

phmatrix create_tridiag_hmatrix (double * a, double * b, double * c, 
			         pccluster rc, pccluster cc)
{
  int n = rc->size;
  phmatrix A = NULL;

  assert (rc->size == cc->size);

  if (rc->son == NULL && cc->son == NULL)
    {
      int i;

      A = new_full_hmatrix (rc, cc);

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

      phmatrix A11 = create_tridiag_hmatrix (a, b, c, rc->son[0], cc->son[0]);
      phmatrix A22 = create_tridiag_hmatrix (a + rc->son[0]->size, 
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

      ref_hmatrix(&A->son[0], A11);
      ref_hmatrix(&A->son[1], A21);
      ref_hmatrix(&A->son[2], A12);
      ref_hmatrix(&A->son[3], A22);

      update_hmatrix (A);
    }

  return A;
}

size_t hmatrix_get_rank (pchmatrix H)
{
  if (H->son)
    return H->son[1]->r->A.cols;
  else
    return H->rc->size;
}

size_t hmatrix_get_m (pchmatrix H)
{
  return H->rc->size;
}

size_t hmatrix_get_n (pchmatrix H)
{
  return H->cc->size;
}
