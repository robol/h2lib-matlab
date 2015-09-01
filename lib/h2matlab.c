#include "hmatlab.h"

extern double ddot(long * n, double *, long *, double *, long *);

#define MATRIX_ELEM(A,i,j,n) A[(j)*(n)+(i)]

ph2matrix create_tridiag_h2matrix (double * a, double * b, double * c, pccluster rc, pccluster cc, pclusterbasis brc, pclusterbasis bcc)
{
  int n = rc->size;
  ph2matrix A = NULL;

  if (brc == NULL)
    {
      brc = create_tridiag_clusterbasis (rc);
      update_tree_clusterbasis (brc);
    }

  if (bcc == NULL)
    {
      bcc = create_tridiag_clusterbasis (cc);
      update_tree_clusterbasis (bcc);
    }

  if (rc->son == NULL && cc->son == NULL)
    {
      int i;

      A = new_full_h2matrix (brc, bcc);

      memset (A->f->a, 0, sizeof (double) * n * n);

      A->f->a[0] = a[0];
      for (i = 0; i < n - 1; i++)
	{
	  MATRIX_ELEM(A->f->a, i+1, i, n) = b[i];
	  MATRIX_ELEM(A->f->a, i + 1, i + 1, n) = a[i+1];
	  MATRIX_ELEM(A->f->a, i, i + 1, n) = c[i];
	}

      update_h2matrix (A);      
    }
  else 
    {
      /* Perform subdivision */
      A = new_super_h2matrix (brc, bcc, 2, 2);

      ph2matrix A11 = create_tridiag_h2matrix (a, b, c, rc->son[0], cc->son[0], brc->son[0], bcc->son[0]);
      ph2matrix A22 = create_tridiag_h2matrix (a + rc->son[0]->size, 
				               b + rc->son[0]->size,
				               c + rc->son[0]->size,
				               rc->son[1], cc->son[1], brc->son[1], bcc->son[1]);

      ph2matrix A12 = new_uniform_h2matrix (brc->son[0], bcc->son[1]);
      ph2matrix A21 = new_uniform_h2matrix (brc->son[1], bcc->son[0]);

      init_zero_amatrix (&A12->u->S, 2, 2);
      init_zero_amatrix (&A21->u->S, 2, 2);

      A12->u->S.a[1] = c[cc->son[0]->size - 1];
      A21->u->S.a[2] = b[rc->son[0]->size - 1];

      ref_h2matrix(&A->son[0], A11);
      ref_h2matrix(&A->son[1], A21);
      ref_h2matrix(&A->son[2], A12);
      ref_h2matrix(&A->son[3], A22);

      update_h2matrix (A);
    }

  return A;
}



pclusterbasis create_tridiag_clusterbasis(pccluster rc)
{
  pclusterbasis clust=NULL;
  int n = rc->size;

  if (rc->sons == 0)
    {
      clust = new_clusterbasis(rc);
      clust->k = 2;
      init_zero_amatrix (&clust->V, rc->size, 2);
      clust->V.a[0] = 1;
      MATRIX_ELEM(clust->V.a, n-1, 1, n) = 1;
    }
  else
    {
      clust=new_clusterbasis(rc);

      clust->son[0] = create_tridiag_clusterbasis (rc->son[0]);
      init_zero_amatrix (&clust->son[0]->E, 2, 2);
      clust->son[0]->E.a[0] = 1;     

      clust->son[1] = create_tridiag_clusterbasis (rc->son[1]);
      init_zero_amatrix (&clust->son[1]->E, 2, 2);
      clust->son[1]->E.a[3] = 1;

      clust->sons = 2;
      clust->k = 2;
  }

  return clust;
}


size_t
h2matrix_get_rank (pch2matrix H)
{
  return 1;
}
