#include <mex.h>
#include <hmatlab.h>

void mexFunction (int nlhs, mxArray * plhs[], 
                  int nrhs, const mxArray * prhs[])
{
    /* We need to have the cluster specified */
    pccluster rowCluster = DESERIALIZE_POINTER (mxGetProperty(prhs[2], 0, "cluster"));
    pccluster colCluster = DESERIALIZE_POINTER (mxGetProperty(prhs[3], 0, "cluster"));
    
    /* TODO: Convert prhs[0] into a matrix A. Simple for the real case, messy
     * for the complex one (since MATLAB storage is not the same of FORTRAN
     * storage... */
#ifdef USE_COMPLEX    
    field * A = 
#else
    field * A = mxGetPr(prhs[1]);
#endif
    
    /* This is implemented inside libhmalab, and performs recursive SVDs
     * to find the correct ranks. */
    phmatrix H = create_hmatrix_from_full (A, rowCluster, colCluster, mxGetN(prhs[1]));
    
    SERIALIZE_POINTER_TO_PROPERTY (prhs[0], "hmatrix", H);
}
        
