#include "mex.h"
#include <avector.h>
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <truncation.h>
#include <hmatlab.h>

static void 
change_sign (phmatrix A)
{
  /* TODO: Use BLAS to improve efficiency */

  if (A->son) 
    {
      int i, j;

      for (i = 0; i < A->csons; i++)
	for (j = 0; j < A->rsons; j++)
	  {
	    change_sign (A->son[i * A->rsons + j]);
	  }
    }

  if (A->f)
    {
      int i, j;

      for (i = 0; i < A->f->rows; i++)
	for (j = 0; j < A->f->cols; j++)
	  A->f->a[i + j * A->f->rows] = - A->f->a[i + j * A->f->rows];
    }

  if (A->r)
    {
      int i, j;

      for (i = 0; i < A->r->A.rows; i++)
	for (j = 0; j < A->r->A.cols; j++)
	  A->r->A.a[i + j * A->r->A.rows] = - A->r->A.a[i + j * A->r->A.rows];
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));
  phmatrix C = clone_hmatrix (A);

  change_sign (C);

  SERIALIZE_POINTER (plhs[0], C);
}
