#include <mex.h>
#include <string.h>
#include <math.h>
static void dqmult(double* a, double* b, double* c)
{
  c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  c[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
  c[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
  c[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
  return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // The rotation angle, in double
    // or single float format
    if(mxGetNumberOfElements(prhs[0])!=4) {
      mexErrMsgTxt("1st input must have length 4"); 
    }
    if(mxGetNumberOfElements(prhs[1])!=4) {
      mexErrMsgTxt("2nd input must have length 4");  
    }
    

    int    ndims;
    mwSize dims[10];
    
    double *d = (double *)mxGetPr(prhs[0]);
    double *d2 = (double *)mxGetPr(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix( 1, 4, mxREAL );
    double *outQuat = mxGetDoubles(plhs[0]);
    
    dqmult(d,d2,outQuat);

    return; 
} 

