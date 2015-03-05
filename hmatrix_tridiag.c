#include "mex.h"
#include <hmatrix.h>
#include <cluster.h>
#include <string.h>

static pcluster create_cluster(int *a, int n)
{
pcluster clust=NULL;
if(b-a<4)
  {
    clust=new_cluster(n,a,0,1);
  }
else
  {
    clust=new_cluster(n,a,2,1);
    cluster * c=clust;
    c->son[0]=create_cluster(a,n/2);
    c->son[1]=create_cluster(a+n/2,n-n/2);
  }
return clust;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
int n=mxGetM(prhs[0]);
int*a=malloc(sizeof(int)*n);
int i;
for(i=0;i<n;i++)
    a[i]=i;
pcluster rc=create_cluster(a,n);
pcluster cc=create_cluster(a,n);
update_cluster(rc);
update_cluster(cc);
phmatrix A=new_hmatrix(rc,cc);

}
