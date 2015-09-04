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

  pclusterbasis cbcopy = clone_clusterbasis (B->cb);
  pclusterbasis rbcopy = clone_clusterbasis (A->rb);

  ph2matrix C = clonestructure_h2matrix (A, rbcopy, cbcopy);
  clear_h2matrix (C);

  pclusteroperator rop = prepare_row_clusteroperator (C->rb, C->cb, tm);
  pclusteroperator cop = prepare_col_clusteroperator (C->rb, C->cb, tm);
  
  addmul_h2matrix (1.0, A, false, B, C, rop, cop, tm, h2lib_eps);

  del_clusteroperator (cop);
  del_clusteroperator (rop);

  SERIALIZE_POINTER (plhs[0], C);
}
