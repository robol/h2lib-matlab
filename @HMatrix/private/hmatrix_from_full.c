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
    int m = mxGetM(prhs[1]);
    int n = mxGetN(prhs[1]);
    int i;
    double * Ar = mxGetPr(prhs[1]);
    double * Ai = mxGetPi(prhs[1]);
    double * A = mxMalloc(n * m * sizeof(double) * 2);

    /* Copy all the contents from the origianal matrix */
    for (i = 0; i < n * m; i++) {
        A[2*i] = Ar[i];
        if (Ai)
            A[2*i+1] = Ai[i];
        else
            A[2*i+1] = 0.0;
    }
#else
    field * A = mxGetPr(prhs[1]);
#endif

    /* This is implemented inside libhmalab, and performs recursive SVDs
     * to find the correct ranks. */
    phmatrix H = create_hmatrix_from_full ((field*) A, rowCluster, colCluster, mxGetN(prhs[1]));

#ifdef USE_COMPLEX
    mxFree(A);
#endif

    SERIALIZE_POINTER_TO_PROPERTY (prhs[0], "hmatrix", H);
}
