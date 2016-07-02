#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include <hmatlab.h>

static void
set_row_to_zero(phmatrix A, int row)
{
    if (A->son == NULL) {
        if (A->f) {
          int i, j;
          for (i = 0; i < A->rc->size; i++) {
              if (A->rc->idx[i] == row) {
                  for (j = 0; j < A->cc->size; j++) {
                      if (A->cc->idx[j] < row) {
                          A->f->a[i + j * A->rc->size] = 0.0;
                      }
                  }
              }
           } 
        }

	if (A->r) {
		int i, j;
	for (i = 0; i < A->rc->size; i++) {
		if (A->rc->idx[i] == row) {
                  if (A->cc->idx[0] < row)
			  for (j = 0; j < A->r->k; j++) {
			      A->r->A.a[i + A->rc->size * j] = 0.0;
			  }
	      }
	   } 
        }
    }

   int i, j;

   for (i = 0; i < A->rsons * A->csons; i++)
         set_row_to_zero(A->son[i], row);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
     int i = *mxGetPr(prhs[1]);
     set_row_to_zero(A, i);
}
