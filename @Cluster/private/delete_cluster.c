#include <mex.h>
#include <hmatlab.h>

void mexFunction (int nlhs, mxArray * plhs[], 
		  int nrhs, const mxArray * prhs[])
{
  pcluster cluster = DESERIALIZE_POINTER (mxGetProperty(prhs[0], 0, "cluster"));
  free (cluster->idx);
  del_cluster (cluster);
}
