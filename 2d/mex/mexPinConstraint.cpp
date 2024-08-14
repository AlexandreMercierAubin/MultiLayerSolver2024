/*==========================================================
 * mexPinConstraint.cpp 
 *
 * To compile type: mex -R2018a mexPinConstraint.cpp
 *========================================================*/

#include "mex.h"
#include <math.h>
#include <cmath>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
        /* check for proper number of arguments */
    if ( nrhs != 3 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nrhs","3 inputs required.");
    }
    if ( nlhs != 3 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nlhs","3 outputs required.");
    }

    double* p = mxGetDoubles(prhs[0]);
    double* pinPos = mxGetDoubles(prhs[1]);
    double* mass = mxGetDoubles(prhs[2]);

    //new vars for output
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double *deltaLambda = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 1, 2, mxREAL );
    double *deltax = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double *constraintC = mxGetDoubles(plhs[2]);

    deltax[0] = pinPos[0]-p[0];
    deltax[1] = pinPos[1]-p[0];
    constraintC[0] = sqrt(l0*l0 + l1*l1);
    deltaLambda[0] = -mass[0]*constraintC[0];
}
