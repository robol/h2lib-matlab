#include "mex.h"
#include <avector.h>
#include <h2matrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include <h2update.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  ph2matrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "h2matrix"));
  ph2matrix B = DESERIALIZE_POINTER (mxGetProperty (prhs[1], 0, "h2matrix"));
  ptruncmode tm = new_releucl_truncmode();
  int i, n = A->cb->t->size;

  pclusterbasis cbcopy = clone_clusterbasis (A->cb);
  pclusterbasis rbcopy = clone_clusterbasis (A->rb);
  ph2matrix C = clone_h2matrix (A, rbcopy, cbcopy);

  pclusteroperator rop = prepare_row_clusteroperator (C->rb, C->cb, tm);
  pclusteroperator cop = prepare_col_clusteroperator (C->rb, C->cb, tm);

  double * ones = malloc (sizeof (double) * A->rb->t->size);
  double * zeros = malloc (sizeof (double) * (A->rb->t->size - 1));

  memset (zeros, 0, sizeof (double) * (A->rb->t->size - 1));

  for (i = 0; i < n; i++) 
    ones[i] = 1.0;
  
  
  ph2matrix I = create_tridiag_h2matrix (ones, zeros, zeros, B->rb->t, B->cb->t, NULL, NULL);
  
  addmul_h2matrix (-1.0, I, false, B, C, rop, cop, tm, h2lib_eps);

  free (ones);
  free (zeros);

  del_h2matrix (I);
  del_clusteroperator (cop);
  del_clusteroperator (rop);

  SERIALIZE_POINTER (plhs[0], C);
}
