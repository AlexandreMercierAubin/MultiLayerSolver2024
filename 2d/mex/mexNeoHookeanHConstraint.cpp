/*==========================================================
 * mexNeoHookeanHConstraint.cpp 
 *
 * To compile type: mex -R2018a mexNeoHookeanHConstraint.cpp
 *========================================================*/

#include "mex.h"
#include <math.h>
#include <cmath>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
        /* check for proper number of arguments */
    if ( nrhs != 9 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nrhs","9 inputs required.");
    }
    if ( nlhs != 6 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nlhs","6 outputs required.");
    }

    // could read like this for cleaner file, but likely a bit slower: double volume = *mxGetDoubles(prhs[6]);

    double *h = mxGetDoubles(prhs[0]);
    double* x = mxGetDoubles(prhs[1]);
    double* DmInv = mxGetDoubles(prhs[2]);
    double* materialMu = mxGetDoubles(prhs[3]);
    double* materialLambda = mxGetDoubles(prhs[4]);
    double area = *mxGetDoubles(prhs[5]);
    double* lambda = mxGetDoubles(prhs[6]);
    double* mass = mxGetDoubles(prhs[7]);
    double gamma = *mxGetDoubles(prhs[8]);

    //new vars for output
    plhs[0] = mxCreateDoubleMatrix( 6, 1, mxREAL );
    double* deltax = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 6, 1, mxREAL );
    double* Minv = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix( 6, 1, mxREAL );
    double* gradC = mxGetDoubles(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double* deltaLambda = mxGetDoubles(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double* alphaTilde = mxGetDoubles(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double* constraintC = mxGetDoubles(plhs[5]);

    double t2 = DmInv[0]*DmInv[3]*x[0];
    double t3 = DmInv[2]*DmInv[1]*x[0];
    double t4 = DmInv[0]*DmInv[3]*x[1];
    double t5 = DmInv[2]*DmInv[1]*x[1];
    double t6 = DmInv[0]*DmInv[3]*x[2];
    double t7 = DmInv[2]*DmInv[1]*x[2];
    double t8 = DmInv[0]*DmInv[3]*x[3];
    double t9 = DmInv[2]*DmInv[1]*x[3];
    double t10 = DmInv[0]*DmInv[3]*x[4];
    double t11 = DmInv[2]*DmInv[1]*x[4];
    double t12 = DmInv[0]*DmInv[3]*x[5];
    double t13 = DmInv[2]*DmInv[1]*x[5];
    double t14 = 1.0/area;
    double t15 = 1.0/(h[0]*h[0]);
    double t16 = 1.0/mass[0];
    double t17 = 1.0/mass[1];
    double t18 = 1.0/mass[2];
    double t19 = 1.0/mass[3];
    double t20 = 1.0/mass[4];
    double t21 = 1.0/mass[5];
    double t22 = 1.0/materialLambda[0];
    double t23 = t2*x[3];
    double t24 = t4*x[2];
    double t25 = t3*x[3];
    double t26 = t5*x[2];
    double t27 = t2*x[5];
    double t28 = t4*x[4];
    double t29 = t3*x[5];
    double t30 = t5*x[4];
    double t31 = t6*x[5];
    double t32 = t8*x[4];
    double t33 = t7*x[5];
    double t34 = t9*x[4];
    double t35 = -t3;
    double t36 = -t5;
    double t37 = -t6;
    double t38 = -t7;
    double t39 = -t8;
    double t40 = -t9;
    double t41 = -t10;
    double t42 = -t12;
    double t49 = t14*t15*t22;
    double t43 = -t23;
    double t44 = -t26;
    double t45 = -t28;
    double t46 = -t29;
    double t47 = -t31;
    double t48 = -t34;
    double t50 = lambda[0]*t49;
    double t52 = t2+t7+t35+t37;
    double t53 = t2+t11+t35+t41;
    double t54 = t4+t9+t36+t39;
    double t55 = t4+t13+t36+t42;
    double t56 = t6+t11+t38+t41;
    double t57 = t8+t13+t40+t42;
    double t51 = -t50;
    double t58 = t52*t52;
    double t59 = t53*t53;
    double t60 = t54*t54;
    double t61 = t55*t55;
    double t62 = t56*t56;
    double t63 = t57*t57;
    double t64 = t21*t58;
    double t65 = t19*t59;
    double t66 = t20*t60;
    double t67 = t17*t62;
    double t68 = t18*t61;
    double t69 = t16*t63;
    double t70 = gamma+t24+t25+t27+t30+t32+t33+t43+t44+t45+t46+t47+t48+t51;
    double t71 = t49+t64+t65+t66+t67+t68+t69;
    double t72 = 1.0/t71;
    deltax[0] = t16*t57*t70*t72;
    deltax[1] = -t17*t56*t70*t72;
    deltax[2] = -t18*t55*t70*t72;
    deltax[3] = t19*t53*t70*t72;
    deltax[4] = t20*t54*t70*t72;
    deltax[5] = -t21*t52*t70*t72;
    Minv[0] = t16;
    Minv[1] = t17;
    Minv[2] = t18;
    Minv[3] = t19;
    Minv[4] = t20;
    Minv[5] = t21;
    gradC[0] = t57;
    gradC[1] = t7+t10-t11+t37;
    gradC[2] = -t4+t5+t12-t13;
    gradC[3] = t53;
    gradC[4] = t54;
    gradC[5] = -t2+t3+t6+t38;
    deltaLambda[0] = t70*t72;
    alphaTilde[0] = t49;
    constraintC[0] = -gamma+t23-t24-t25+t26-t27+t28+t29-t30+t31-t32-t33+t34;

}
