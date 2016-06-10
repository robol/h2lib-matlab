#include "hmatlab.h"
#include <complex.h>
#include <string.h>
#include <cluster.h>
#define MATRIX_ELEM(A,i,j,n) A[(j)*(n)+(i)]

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

/**
 * @brief Create a hierarichal matrix representing a tridiagonal
 * matrix.
 *
 * @param a The diagonal entries of the matrix.
 * @param b The subdiagonal entries of the matrix.
 * @param c The superdiagonal entries of the matrix.
 */
phmatrix create_real_tridiag_hmatrix (double * a, double * b, double * c,
				 	        	      pccluster rc, pccluster cc)
{
  field *a_a = malloc (rc->size * sizeof (field));
  field *b_b = malloc (rc->size * sizeof (field));
  field *c_c = malloc (rc->size * sizeof (field));
  int i;

  for (i = 0; i < rc->size; i++) {
    a_a[i] = a[i];
    b_b[i] = b[i];
    c_c[i] = c[i];
  }

 phmatrix T = create_tridiag_hmatrix (a_a, b_b, c_c, rc, cc);

 free (a_a);
 free (b_b);
 free (c_c);

 return T;
}

#ifndef USE_COMPLEX
extern void
dgesvd_ (const char *, const char *, LAPACK_INT*, LAPACK_INT*,
	 field *, LAPACK_INT*, double *, field *, LAPACK_INT*,
	 field *, LAPACK_INT*, field*, LAPACK_INT*, LAPACK_INT*);
#else
extern void
zgesvd_ (const char *, const char *, LAPACK_INT*, LAPACK_INT*,
	 field*, LAPACK_INT*, double *, field *, LAPACK_INT*,
	 field*, LAPACK_INT*, field*, LAPACK_INT*, double*, LAPACK_INT*);
#endif

phmatrix
constructCompressedRkMatrix(field * X, pccluster rc, pccluster cc, double eps)
{
  int rank = 0, i, j;
  LAPACK_INT m = rc->size;
  LAPACK_INT n = cc->size;
  LAPACK_INT mn = (m < n) ? m : n;
  LAPACK_INT MN = (m > n) ? m : n;
  phmatrix H = NULL;
  LAPACK_INT info = 0, lwork = 10 * MN;
  field * work = malloc (sizeof (field) * lwork);
  field * U = malloc (sizeof (field) * m * m);
  field * Vt = malloc (sizeof (field) * n * n);

  double * singular_values = malloc (sizeof (double) * m);

  /* Make a workspace query */
#ifndef USE_COMPLEX
  dgesvd_("S", "S", &m, &n, X, &m, singular_values,
	  U, &m, Vt, &n, work, &lwork, &info);
#else
  double * rwork = malloc(sizeof(double) * 5 * mn);
  zgesvd_("S", "S", &m, &n, X, &m, singular_values,
	  U, &m, Vt, &n, work, &lwork, rwork, &info);
  free(rwork);
#endif

  /* Find out the rank */
  while (singular_values[rank] > eps && rank < mn) {
    rank++;
  }

  /* Construct the matrix */
  H = new_rk_hmatrix (rc, cc, rank);

  /* Fill in the basis */
  memcpy (H->r->A.a, U,  sizeof(field) * rank * m);
  for (i = 0; i < rank; i++)
    for (j = 0; j < n; j++)
#ifdef USE_COMPLEX
      MATRIX_ELEM(H->r->B.a, j, i, n) = conj(MATRIX_ELEM(Vt, i, j, n)) * singular_values[i];
#else
      MATRIX_ELEM(H->r->B.a, j, i, n) = conj(MATRIX_ELEM(Vt, i, j, n)) * singular_values[i];
#endif

  free(U);
  free(Vt);
  free(singular_values);

  return H;
}

phmatrix create_hmatrix_from_full (field * a, pccluster rc, pccluster cc, int lda)
{
  phmatrix A = NULL;
  int n = rc->size, i, j;

  /* Check if we are in the base case */
  if (rc->son == NULL && cc->son == NULL) {
    A = new_full_hmatrix (rc, cc);

    /* Copy the content inside the full matrix */
    for (i = 0; i < rc->size; i++)
      for (j = 0; j < cc->size; j++)
	MATRIX_ELEM(A->f->a, i, j, n) = MATRIX_ELEM(a, rc->idx[i], cc->idx[j], lda);
  }
  else {
    A = new_super_hmatrix (rc, cc, 2, 2);

    phmatrix A11 = create_hmatrix_from_full (a, rc->son[0], cc->son[0], lda);
    phmatrix A22 = create_hmatrix_from_full (a, rc->son[1], cc->son[1], lda);

    /* We need to determine the rank of the offdiagonal blocks. */
    field * X = malloc (sizeof(field) * rc->son[1]->size * cc->son[0]->size);
    for (i = 0; i < rc->son[1]->size; i++)
      for (j = 0; j < cc->son[0]->size; j++)
      MATRIX_ELEM(X, i, j, rc->son[1]->size) = MATRIX_ELEM(a, rc->son[1]->idx[i], cc->son[0]->idx[j], lda);

    /* Determine the rank of the matrix */
    phmatrix A21 = constructCompressedRkMatrix(X, rc->son[1], cc->son[0], h2lib_eps);

    /* Try to rellocate X if possible */
    field * newX = realloc(X, sizeof(field) * rc->son[0]->size * cc->son[1]->size);
    if (newX != NULL)
    X = newX;

    for (i = 0; i < rc->son[0]->size; i++) {
      for (j = 0; j < cc->son[1]->size; j++) {
	MATRIX_ELEM(X, i, j, rc->son[0]->size) = MATRIX_ELEM(a, rc->son[0]->idx[i], cc->son[1]->idx[j], lda);
      }
    }
    phmatrix A12 = constructCompressedRkMatrix(X, rc->son[0], cc->son[1], h2lib_eps);
    free(X);

    ref_hmatrix(&A->son[0], A11);
    ref_hmatrix(&A->son[1], A21);
    ref_hmatrix(&A->son[2], A12);
    ref_hmatrix(&A->son[3], A22);

    update_hmatrix (A);
  }

  return A;
}


phmatrix create_tridiag_hmatrix (field * a, field * b, field * c,
			         			 pccluster rc, pccluster cc)
{
  int n = rc->size;
  phmatrix A = NULL;

  assert (rc->size == cc->size);

  if (rc->son == NULL && cc->son == NULL)
    {
      int i;

      A = new_full_hmatrix (rc, cc);

      memset (A->f->a, 0, sizeof (field) * n * n);

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

      memset (A12->r->A.a, 0, sizeof (field) * rc->son[0]->size);
      memset (A12->r->B.a, 0, sizeof (field) * cc->son[1]->size);
      memset (A21->r->A.a, 0, sizeof (field) * rc->son[1]->size);
      memset (A21->r->B.a, 0, sizeof (field) * cc->son[0]->size);

      /* Fill in the low rank parts */
      A12->r->A.a[rc->son[0]->size - 1] = c[cc->son[0]->size - 1];
      A12->r->B.a[0] = 1.0;
      A21->r->A.a[0] = b[rc->son[0]->size - 1];
      A21->r->B.a[cc->son[0]->size - 1] = 1.0;

      ref_hmatrix(&A->son[0], A11);
      ref_hmatrix(&A->son[1], A21);
      ref_hmatrix(&A->son[2], A12);
      ref_hmatrix(&A->son[3], A22);

      update_hmatrix (A);
    }

  return A;
}

pcluster shift_cluster(pccluster original_cluster, int shift) 

{ 
int i;

pcluster newcluster = new_cluster(original_cluster->size, original_cluster->idx + shift, original_cluster->sons, original_cluster->dim);

for ( i = 0; i < original_cluster-> sons; i++)
		newcluster->son[i]=shift_cluster(original_cluster->son[i], shift);

return newcluster;

}


phmatrix create_band_hmatrix (double * a, double * b, double * c, int p, int q,
			         pccluster rc, pccluster cc)
{
  int n = rc->size;
  int i,j;
  phmatrix A = NULL;

  assert (rc->size == cc->size);

  if (rc->son == NULL && cc->son == NULL)
    {


      A = new_full_hmatrix (rc, cc);

      memset (A->f->a, 0, sizeof (double) * n * n);

      for (i = 0; i < n ; i++)
	{
          for (j = 0; j < i; j++)
            {
              if (i-j <= p)
	        {
                  MATRIX_ELEM(A->f->a, i, j, n) = b[p*j+i-j-1];
                }
            }
	  MATRIX_ELEM(A->f->a, i, i, n) = a[i];
          for (j = i+1; j < n; j++)
	    {
            if (j-i <= q)
	        {
                  MATRIX_ELEM(A->f->a, i, j, n) = c[q*i+j-i-1];
                }
            }
	}

      update_hmatrix (A);
    }
  else
    {
      /* Perform subdivision */
      A = new_super_hmatrix (rc, cc, 2, 2);

      phmatrix A11 = create_band_hmatrix (a, b, c, p, q, rc->son[0], cc->son[0]);
      phmatrix A22 = create_band_hmatrix (a + rc->son[0]->size,
				             b + rc->son[0]->size * p,
				             c + rc->son[0]->size * q,
				             p, q, rc->son[1], cc->son[1]);

      phmatrix A12 = new_rk_hmatrix (rc->son[0], cc->son[1], q);
      phmatrix A21 = new_rk_hmatrix (rc->son[1], cc->son[0], p);

      memset (A12->r->A.a, 0, sizeof (field) * rc->son[0]->size *q);
      memset (A12->r->B.a, 0, sizeof (field) * cc->son[1]->size *q);
      memset (A21->r->A.a, 0, sizeof (field) * rc->son[1]->size *p);
      memset (A21->r->B.a, 0, sizeof (field) * cc->son[0]->size *p);

      /* Fill in the low rank parts */
      for (i=0; i<q; i++)
        {
          A12->r->A.a[(i+1) * (rc->son[0]->size - 1)] = 1.0;
          for (j=0; j<q-i; j++)
            {
              A12->r->B.a[i * cc->son[1]->size + j] = c[(rc->son[0]->size - i - 1)* q + j + i];
            }
        }
      for (i=0; i<p; i++)
        {
          A21->r->A.a[i* (rc->son[1]->size + 1)] = 1.0;
          for (j=cc->son[0]->size - p + i; j < cc->son[0]->size; j++)
            {
              A21->r->B.a[i * cc->son[0]->size + j] = b[p * j + cc->son[0]->size - j - 1 + i];
            }
        }
      ref_hmatrix(&A->son[0], A11);
      ref_hmatrix(&A->son[1], A21);
      ref_hmatrix(&A->son[2], A12);
      ref_hmatrix(&A->son[3], A22);

      update_hmatrix (A);
    }

  return A;
}

phmatrix create_generators_hmatrix (double * d, double * U, double * V,  double * W, double * Z,
			         int ksub,  int ksup, pccluster rc, pccluster cc)
{
  int n = rc->size;
  phmatrix A = NULL;
  int i, j, t;

  assert (rc->size == cc->size);

  if (rc->son == NULL && cc->son == NULL)
    {
      A = new_full_hmatrix (rc, cc);

      memset (A->f->a, 0, sizeof (field) * n * n);

      for (i = 0; i < n; i++)
	{
          for (j=0; j < i; j++)
            {
              for (t=0; t < ksub; t++)
                {
                  MATRIX_ELEM(A->f->a, i, j, n) += U[i * ksub + t] * V[j * ksub + t];
                }
            }
          for (j=i+1; j < n; j++)
            {
              for (t=0; t < ksup; t++)
                {
                  MATRIX_ELEM(A->f->a, i, j, n) += W[i * ksup + t] * Z[j * ksup + t];
                }
            }
          MATRIX_ELEM(A->f->a, i, i, n) = d[i];
	}

      update_hmatrix (A);
    }
  else
    {
      /* Perform subdivision */
      A = new_super_hmatrix (rc, cc, 2, 2);

      phmatrix A11 = create_generators_hmatrix (d,U,V,W,Z,ksub,ksup, rc->son[0], cc->son[0]);
      phmatrix A22 = create_generators_hmatrix (d + rc->son[0]->size,
				             U + rc->son[0]->size * ksub,
				             V + cc->son[0]->size * ksub,
				             W + rc->son[0]->size * ksup,
				             Z + cc->son[0]->size * ksup,
				             ksub, ksup, rc->son[1], cc->son[1]);

      phmatrix A12 = new_rk_hmatrix (rc->son[0], cc->son[1], ksup);
      phmatrix A21 = new_rk_hmatrix (rc->son[1], cc->son[0], ksub);
      memset (A12->r->A.a, 0, sizeof (field) * rc->son[0]->size * ksup);
      memset (A12->r->B.a, 0, sizeof (field) * cc->son[1]->size * ksup);
      memset (A21->r->A.a, 0, sizeof (field) * rc->son[1]->size * ksub);
      memset (A21->r->B.a, 0, sizeof (field) * cc->son[0]->size * ksub);

      /* Fill in the low rank parts */
      for (i=0; i < rc->son[0]->size; i++)
        {
          for(j=0; j < ksup; j++)
            {
              A12->r->A.a[j * rc->son[0]->size + i] = W[i * ksup +j];
            }
        }
      for (i=0; i < cc->son[1]->size; i++)
        {
          for(j=0; j < ksup; j++)
            {
              A12->r->B.a[j * cc->son[1]->size + i] = Z[cc->son[0]->size * ksup + i * ksup + j];
            }
        }
      for (i=0; i < rc->son[1]->size; i++)
        {
          for(j=0; j < ksub; j++)
            {
              A21->r->A.a[j * rc->son[1]->size + i] = U[rc->son[0]->size * ksub + i * ksub + j];
            }
        }
      for (i=0; i < cc->son[0]->size; i++)
        {
          for(j=0; j < ksub; j++)
            {
              A21->r->B.a[j * cc->son[0]->size + i] = V[i * ksub +j];
            }
        }
      ref_hmatrix(&A->son[0], A11);
      ref_hmatrix(&A->son[1], A21);
      ref_hmatrix(&A->son[2], A12);
      ref_hmatrix(&A->son[3], A22);

      update_hmatrix (A);
    }

  return A;
}

field hmatrix_trace(pchmatrix H)
{
if (H->f){
	return trace_amatrix(H->f);
	}
if(H->son){
	return hmatrix_trace(H->son[0])+hmatrix_trace(H->son[3]);
	}
}

field
trace_amatrix(pamatrix a)
{
  LAPACK_INT      rows = a->rows;
  LAPACK_INT lda = a->ld;
  field trace = 0.0;
  uint      j;

  for (j = 0; j < rows; j++) {
    trace+=a->a[j + j * lda];
  }
return trace;
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

void scale_hmatrix (field r, phmatrix H)
{
  if (H->f) {
    scale_amatrix (r, H->f);
  }

  if (H->r) {
    scale_amatrix (r, &H->r->A);
  }

  if (H->son) {
    int i, j;

    for (i = 0; i < H->rsons * H->csons; i++)
      scale_hmatrix (r, H->son[i]);
  }
}

void shift_hmatrix (field r, phmatrix H)
{
  if (H->f) {
    int i;

    for (i = 0; i < H->rc->size; i++)
      setentry_amatrix (H->f, i, i, getentry_amatrix (H->f, i, i) + r);
  }

  if (H->son) {
    int i;

    for (i = 0; i < H->rsons * H->csons; i++)
      shift_hmatrix (r, H->son[i]);
  }
}
