/*==========================================================
 * mexSTVKConstraintEigen.cpp 
 *
 * To compile type: mex -R2018a mexStVKVoigtConstraintEigen.cpp  -I"../../lib/gptoolbox/mex/external/libigl/external/eigen" -I"../../lib/gptoolbox/mex/external/libigl/include"
 *========================================================*/

#include "mex.h"
#include <math.h>
#include <cmath>
#include <Eigen/Dense>

#define NORMALIZE_GRAD false

using namespace Eigen;

inline void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
        /* check for proper number of arguments */
    if ( nrhs != 5 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nrhs","5 inputs required.");
    }
    if ( nlhs != 4 && nlhs != 2) {
        mexErrMsgIdAndTxt("ARP:constraint:nlhs","4 or 2 outputs required.");
    }

    // could read like this for cleaner file, but likely a bit slower: double volume = *mxGetDoubles(prhs[6]);

    double *x = mxGetDoubles(prhs[0]);
    double *DmInv = mxGetDoubles(prhs[1]);
    double *alphaTilde_in = mxGetDoubles(prhs[2]);
    double *lambda_in = mxGetDoubles(prhs[3]);
    double *massInv_in = mxGetDoubles(prhs[4]);
   
    //new vars for output
    plhs[0] = mxCreateDoubleMatrix( 3, 1, mxREAL );
    double *constraintC = mxGetDoubles(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix( 18, 1, mxREAL );
    double *gradC = mxGetDoubles(plhs[1]);

    double *deltaLambda_out;
    double *deltax_out;
    if (nlhs != 2){
        plhs[2] = mxCreateDoubleMatrix( 3, 1, mxREAL );
        deltaLambda_out = mxGetDoubles(plhs[2]);
    
        plhs[3] = mxCreateDoubleMatrix( 6, 1, mxREAL );
        deltax_out = mxGetDoubles(plhs[3]);
    }
    

    double t2 = DmInv[0]+DmInv[1];
    double t3 = DmInv[2]+DmInv[3];
    double t4 = -x[4];
    double t5 = -x[5];
    double t6 = t4+x[0];
    double t7 = t5+x[1];
    double t8 = t4+x[2];
    double t9 = t5+x[3];
    double t10 = DmInv[0]*t6;
    double t11 = DmInv[2]*t6;
    double t12 = DmInv[0]*t7;
    double t13 = DmInv[2]*t7;
    double t14 = DmInv[1]*t8;
    double t15 = DmInv[3]*t8;
    double t16 = DmInv[1]*t9;
    double t17 = DmInv[3]*t9;
    double t18 = t10+t14;
    double t19 = t11+t15;
    double t20 = t12+t16;
    double t21 = t13+t17;
    gradC[0] = DmInv[0]*t18;
    gradC[1] = DmInv[0]*t20;
    gradC[2] = DmInv[1]*t18;
    gradC[3] = DmInv[1]*t20;
    gradC[4] = -t2*t18;
    gradC[5] = -t2*t20;
    gradC[6] = DmInv[2]*t19;
    gradC[7] = DmInv[2]*t21;
    gradC[8] = DmInv[3]*t19;
    gradC[9] = DmInv[3]*t21;
    gradC[10] = -t3*t19;
    gradC[11] = -t3*t21;
    gradC[12] = (DmInv[0]*t19)/2.0+(DmInv[2]*t18)/2.0;
    gradC[13] = (DmInv[0]*t21)/2.0+(DmInv[2]*t20)/2.0;
    gradC[14] = (DmInv[1]*t19)/2.0+(DmInv[3]*t18)/2.0;
    gradC[15] = (DmInv[1]*t21)/2.0+(DmInv[3]*t20)/2.0;
    gradC[16] = t2*t19*(-1.0/2.0)-(t3*t18)/2.0;
    gradC[17] = t2*t21*(-1.0/2.0)-(t3*t20)/2.0;
    constraintC[0] = (t18*t18)/2.0+(t20*t20)/2.0-1.0/2.0;
    constraintC[1] = (t19*t19)/2.0+(t21*t21)/2.0-1.0/2.0;
    constraintC[2] = (t18*t19)/2.0+(t20*t21)/2.0;

    if(NORMALIZE_GRAD){
        for(size_t i = 0; i<8 ; ++i){
            double normVert = sqrt(gradC[i*2]*gradC[i*2] + gradC[i*2+1]*gradC[i*2+1]);
            gradC[i*2] /= normVert;
            gradC[i*2+1] /= normVert;
        }
    }

    if(nlhs == 2){return;}

    MatrixXd alphaTilde = Map<MatrixXd>( alphaTilde_in, 3, 3 );
    MatrixXd lambda = Map<MatrixXd>( lambda_in, 3, 1 );
    MatrixXd massInv = Map<MatrixXd>( massInv_in, 6, 1 );
    MatrixXd Minv = massInv.asDiagonal();
    MatrixXd gradCmat = Map<MatrixXd>( gradC, 6, 3 );
    MatrixXd C = Map<MatrixXd>( constraintC, 3, 1 );

    MatrixXd denominator = gradCmat.transpose()*Minv*gradCmat + alphaTilde;
    MatrixXd numerator = -C - alphaTilde*lambda;
    MatrixXd deltaLambda = denominator.llt().solve(numerator);
    MatrixXd deltax = Minv*gradCmat*deltaLambda;

    for(size_t i = 0; i<3; ++i){
        deltax_out[i*2] = deltax(i*2);
        deltax_out[i*2+1] = deltax(i*2+1);
        deltaLambda_out[i] = deltaLambda(i);
    }
}
