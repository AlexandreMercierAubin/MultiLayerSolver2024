/*==========================================================
 * mexRigidBoundaryConstraint.cpp 
 *
 * To compile type: mex -R2018a mexRigidBoundaryConstraint.cpp
 *========================================================*/

#include "mex.h"
#include <math.h>
#include <cmath>

double cross2D(double v11, double v12, double v21, double v22){
    return v11*v22 - v12*v21;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
        /* check for proper number of arguments */
    if ( nrhs != 12 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nrhs","12 inputs required.");
    }
    if ( nlhs != 4 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nlhs","4 outputs required.");
    }

    double* h = mxGetDoubles(prhs[0]);
    double* pElastic = mxGetDoubles(prhs[1]);
    double* pRigid = mxGetDoubles(prhs[2]);
    double* r = mxGetDoubles(prhs[3]);
    double* massElastic = mxGetDoubles(prhs[4]);
    double* lambda = mxGetDoubles(prhs[5]);
    double* alpha = mxGetDoubles(prhs[6]);
    double* rigidBodyMass = mxGetDoubles(prhs[7]);
    double* iterateCOM = mxGetDoubles(prhs[8]);
    double* iterateRotation = mxGetDoubles(prhs[9]);
    double* iterateInertia = mxGetDoubles(prhs[10]);
    double* nb = mxGetDoubles(prhs[11]);

    //new vars for output
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double *deltaLambda = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 1, 2, mxREAL );
    double *deltaxElastic = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix( 1, 2, mxREAL );
    double *rigidCOM = mxGetDoubles(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double *rigidRotation = mxGetDoubles(plhs[3]);

    double minv = 1.0/massElastic[0];
    double d0 = pElastic[0]-pRigid[0];
    double d1 = pElastic[1]-pRigid[1];
    double constraintC = sqrt(d0*d0 + d1*d1);

    if(abs(constraintC) <= 0.000000001){
        deltaLambda[0] = 0;
        rigidRotation[0] = iterateRotation[0];
        deltaxElastic[0] = 0;
        deltaxElastic[1] = 0;
        rigidCOM[0] = iterateCOM[0];
        rigidCOM[1] = iterateCOM[1];
        return;
    }

    double Iinv = 1/iterateInertia[0];
    double invMassRigid = 1/rigidBodyMass[0];
    double invMassRigidScaled = invMassRigid/nb[0];

    double n0 = d0/constraintC;
    double n1 = d1/constraintC;

    double crossrn= cross2D(r[0],r[1],n0,n1);
    double rigidBodyKineticAtPoint = Iinv*crossrn*crossrn;
    double w = invMassRigidScaled + rigidBodyKineticAtPoint;
    double minvSum = minv+w;
    double h2 = h[0]*h[0];
    double alphaTilde = alpha[0]/h2;

    double numerator = -constraintC - alphaTilde*lambda[0];
    double denumerator = minvSum + alphaTilde;

    deltaLambda[0] = numerator / denumerator;

    double p0 = deltaLambda[0]*n0;
    double p1 = deltaLambda[0]*n1;
    
    deltaxElastic[0] = minv*p0;
    deltaxElastic[1] = minv*p1;

    rigidCOM[0] = iterateCOM[0] - (p0*invMassRigid);
    rigidCOM[1] = iterateCOM[1] - (p1*invMassRigid);

    double deltaOmega = 0.5*Iinv*cross2D(r[0],r[1],p0,p1);
    rigidRotation[0] = iterateRotation[0] - deltaOmega;
}
