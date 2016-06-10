#include "mex.h"
#include <hmatrix.h>
#include <hmatlab.h>
void replace_clusters(phmatrix A, pccluster rc, pccluster cc) {
  int i, j;

  A->rc = rc;
  A->cc = cc;

  for (i = 0; i < A->rsons; i++) {
    for (j = 0; j < A->csons; j++) {
      replace_clusters(A->son[(j)*A->rsons + i], rc->son[i], cc->son[j]);
    }
}
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  phmatrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "hmatrix"));

  /* Get the leading block of A, and the associated cluster trees. */
  phmatrix Al = NULL;
  ref_hmatrix(&Al, clone_hmatrix(A->son[3]));
 
   
 pcluster lrc = shift_cluster(Al->rc, -Al->rc->idx[0]); 
 update_cluster(lrc); 
 pcluster lcc = shift_cluster(Al->cc, -Al->cc->idx[0]);
update_cluster(lcc); 
 replace_clusters(Al,lrc,lcc);
  
  
  /* Copy the relevant data in the output */
  SERIALIZE_POINTER_TO_PROPERTY(prhs[1], "hmatrix", Al);
  SERIALIZE_POINTER_TO_PROPERTY(prhs[2], "cluster", lrc);
  SERIALIZE_POINTER_TO_PROPERTY(prhs[3], "cluster", lcc);
} 


