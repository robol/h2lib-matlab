#include "mex.h"
#include <avector.h>
#include <h2matrix.h>
#include <block.h>
#include <cluster.h>
#include <clusterbasis.h>
#include <string.h>
#include <truncation.h>
#include <hmatlab.h>
#include <h2arith.h>


static bool massei_robol_admissible_cluster (pcluster rc, pcluster cc, void * data)
{
  if (rc->idx[0] != cc->idx[0])
    return true;
  else
    return false;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     int n = mxGetN(prhs[1]);
     int i;
     int dims[] = { 1, n };

     ph2matrix A = DESERIALIZE_POINTER (mxGetProperty (prhs[0], 0, "h2matrix"));
     
     pclusterbasis Arb = clone_clusterbasis (A->rb);
     pclusterbasis Acb = clone_clusterbasis (A->cb);
     ph2matrix Alr = clone_h2matrix(A, Arb, Acb);

     ptruncmode tm = new_releucl_truncmode();

     pclusterbasis Lrb = build_from_cluster_clusterbasis (A->rb->t);
     pclusterbasis Lcb = build_from_cluster_clusterbasis (A->cb->t);

     pcblock Lb = build_strict_block ((pcluster) Lrb->t, (pcluster) Lcb->t, 
				      NULL, &massei_robol_admissible_cluster);
     ph2matrix L = build_from_block_lower_h2matrix (Lb, Lrb, Lcb);

     pclusterbasis Rrb = build_from_cluster_clusterbasis (A->rb->t);
     pclusterbasis Rcb = build_from_cluster_clusterbasis (A->cb->t);

     pblock Rb = build_strict_block ((pcluster) Rrb->t, (pcluster) Rcb->t, 
				     NULL, &massei_robol_admissible_cluster);
     ph2matrix R = build_from_block_upper_h2matrix (Rb, Rrb, Rcb);

     pclusteroperator Arop = prepare_row_clusteroperator (Arb, Acb, tm);
     pclusteroperator Acop = prepare_col_clusteroperator (Arb, Acb, tm);

     pclusteroperator Lrop = prepare_row_clusteroperator (L->rb, L->cb, tm);
     pclusteroperator Lcop = prepare_col_clusteroperator (L->rb, L->cb, tm);

     pclusteroperator Rrop = prepare_row_clusteroperator (R->rb, R->cb, tm);
     pclusteroperator Rcop = prepare_col_clusteroperator (R->rb, R->cb, tm);
    
     lrdecomp_h2matrix(Alr, Arop, Acop, L, Lrop, Lcop, R, Rrop, Rcop, tm, h2lib_eps);

     plhs[0] = mxCreateNumericArray (2, dims, mxINT64_CLASS, mxREAL);
     long * pointers = (long*) mxGetData(plhs[0]);

     del_h2matrix (Alr);
     del_clusteroperator (Arop);
     del_clusteroperator (Acop);

     for (i = 0; i < n; i++)
       {
	 pclusterbasis temprb = clone_clusterbasis (A->rb);
	 pclusterbasis tempcb = clone_clusterbasis (A->cb);
	 ph2matrix temp = clonestructure_h2matrix (A, temprb, tempcb);
         pclusteroperator temprop = prepare_row_clusteroperator (temp->rb, temp->cb, tm);
         pclusteroperator tempcop = prepare_col_clusteroperator (temp->rb, temp->cb, tm);

	 clear_h2matrix (temp);

         ph2matrix B = DESERIALIZE_POINTER (mxGetProperty (prhs[1], i, "h2matrix"));
         pclusterbasis Brb = clone_clusterbasis (B->rb);
         pclusterbasis Bcb = clone_clusterbasis (B->cb);
         ph2matrix Blr = clone_h2matrix(B, Brb, Bcb);

         pclusteroperator Brop = prepare_row_clusteroperator (Blr->rb, Blr->cb, tm);
         pclusteroperator Bcop = prepare_col_clusteroperator (Blr->rb, Blr->cb, tm);

         lowersolve_h2matrix (true, false, L, false, Blr, Brop, Bcop, 
	    temp, temprop, tempcop, tm, h2lib_eps);
         uppersolve_h2matrix (false, false, R, false, 
	    temp, temprop, tempcop, tm, h2lib_eps);

	 pointers[i] = (long) temp;

         del_clusteroperator (Brop);
         del_clusteroperator (Bcop);
	 del_clusteroperator (temprop);
	 del_clusteroperator (tempcop);
         del_h2matrix (Blr);
       }

     del_h2matrix(L);
     del_h2matrix(R);
}
