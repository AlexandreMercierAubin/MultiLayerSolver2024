#include "mex.h"
#include "blas.h"
#include <algorithm>
#include <cmath>
#include "solveTools.h"

//only if using the Eigen version
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <Eigen/Core>
#include <Eigen/Dense>
//mex -R2018a mexSolbe3by3.cpp  -I"../../lib/gptoolbox/mex/external/libigl/external/eigen" -I"../../lib/gptoolbox/mex/external/libigl/include"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    if ( nrhs != 2 ) {
        mexErrMsgIdAndTxt("ARP:mexSolve3by3:nrhs","2 inputs required.");
    }
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 9 ) {
        mexErrMsgIdAndTxt("ARP:mexSolve3by3:size","A must be 3x3 double matrix");
    }
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexSolve3by3:size","b must be size 3.");
    }

//     double *A = (double *)mxGetPr(prhs[0]);
//     double *b = (double *)mxGetPr(prhs[0]);
//     plhs[0] = mxCreateDoubleMatrix( 3, 1, mxREAL );
//     double *x = mxGetDoubles(plhs[0]);
//     solve3by3(A,b,x);

    Eigen::MatrixXd A;
    Eigen::VectorXd b, x;
	
	igl::matlab::parse_rhs_double(prhs+0,A);
	igl::matlab::parse_rhs_double(prhs+1,b);
    
    x = A.llt().solve(b);
	
	igl::matlab::prepare_lhs_double(x,plhs+0);

    return;
}
