#include "mex.h"
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>
#include <assert.h>
#include <hmatlab.h>
#include <complex.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  pcluster rc = DESERIALIZE_POINTER (mxGetProperty (prhs[1], 0, "cluster"));
  pcluster cc = DESERIALIZE_POINTER (mxGetProperty (prhs[2], 0, "cluster"));
  
  field *a, *b, *c;
  int i, n;
  
  n = mxGetM(prhs[3]) > mxGetN(prhs[3]) ? mxGetM(prhs[3]) : mxGetN(prhs[3]);
  
  double *r_a = mxGetPr(prhs[3]), *i_a = mxGetPi(prhs[3]);
  double *r_b = mxGetPr(prhs[4]), *i_b = mxGetPi(prhs[4]);
  double *r_c = mxGetPr(prhs[5]), *i_c = mxGetPi(prhs[5]);    
  
  a = malloc (sizeof (field) * n);
  b = malloc (sizeof (field) * n);
  c = malloc (sizeof (field) * n);
  
  for (i = 0; i < n; i++)
  {
  	a[i] = r_a[i] + I * (i_a ? i_a[i] : 0);
  	b[i] = r_b[i] + I * (i_b ? i_b[i] : 0);
  	c[i] = r_c[i] + I * (i_c ? i_c[i] : 0);
  }

  phmatrix A = create_tridiag_hmatrix (a, b, c, rc, cc);
				  
  free (a);
  free (b);
  free (c);

  SERIALIZE_POINTER_TO_PROPERTY (prhs[0], "hmatrix", A);
}
