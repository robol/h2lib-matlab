#include <mex.h>
#include <cluster.h>
#include <hmatrix.h>
#include <hmatlab.h>

void mexFunction(int nlhs, mxArray* plhs[],
		 int nrhs, const mxArray* prhs[])
{
  pccluster rowCluster = DESERIALIZE_POINTER (mxGetProperty(prhs[1], 0, "cluster"));
  pccluster colCluster = DESERIALIZE_POINTER (mxGetProperty(prhs[2], 0, "cluster"));

  phmatrix H = new_hmatrix (rowCluster, colCluster);
  
  SERIALIZE_POINTER_TO_PROPERTY (prhs[0], "hmatrix", H);
}
