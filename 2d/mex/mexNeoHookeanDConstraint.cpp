/*==========================================================
 * mexNeoHookeanDConstraint.cpp 
 *
 * To compile type: mex -R2018a mexNeoHookeanDConstraint.cpp
 *========================================================*/

#include "mex.h"
#include <math.h>
#include <cmath>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
        /* check for proper number of arguments */
    if ( nrhs != 8 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nrhs","8 inputs required.");
    }
    if ( nlhs != 6 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nlhs","6 outputs required.");
    }

    // could read like this for cleaner file, but likely a bit slower: double volume = *mxGetDoubles(prhs[6]);

    double* h = mxGetDoubles(prhs[0]);
    double* x = mxGetDoubles(prhs[1]);
    double* DmInv = mxGetDoubles(prhs[2]);
    double* materialMu = mxGetDoubles(prhs[3]);
    double* materialLambda = mxGetDoubles(prhs[4]);
    double area = *mxGetDoubles(prhs[5]);
    double* lambda = mxGetDoubles(prhs[6]);
    double* mass = mxGetDoubles(prhs[7]);

    //new vars for output
    plhs[0] = mxCreateDoubleMatrix( 6, 1, mxREAL );
    double *deltax = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 6, 1, mxREAL );
    double *Minv = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix( 6, 1, mxREAL );
    double *gradC = mxGetDoubles(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double *deltaLambda = mxGetDoubles(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double *alphaTilde = mxGetDoubles(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double *constraintC = mxGetDoubles(plhs[5]);
    
    double t2 = DmInv[0]+DmInv[1];
    double t3 = DmInv[2]+DmInv[3];
    double t4 = 1.0/area;
    double t5 = 1.0/(h[0]*h[0]);
    double t6 = 1.0/mass[0];
    double t7 = 1.0/mass[1];
    double t8 = 1.0/mass[2];
    double t9 = 1.0/mass[3];
    double t10 = 1.0/mass[4];
    double t11 = 1.0/mass[5];
    double t12 = 1.0/materialMu[0];
    double t13 = -x[4];
    double t14 = -x[5];
    double t15 = t13+x[0];
    double t16 = t14+x[1];
    double t17 = t13+x[2];
    double t18 = t14+x[3];
    double t27 = t4*t5*t12;
    double t19 = DmInv[0]*t15;
    double t20 = DmInv[2]*t15;
    double t21 = DmInv[0]*t16;
    double t22 = DmInv[2]*t16;
    double t23 = DmInv[1]*t17;
    double t24 = DmInv[3]*t17;
    double t25 = DmInv[1]*t18;
    double t26 = DmInv[3]*t18;
    double t28 = lambda[0]*t27;
    double t29 = t19+t23;
    double t30 = t20+t24;
    double t31 = t21+t25;
    double t32 = t22+t26;
    double t33 = t29*t29;
    double t34 = t30*t30;
    double t35 = t31*t31;
    double t36 = t32*t32;
    double t37 = DmInv[0]*t29*2.0;
    double t38 = DmInv[2]*t30*2.0;
    double t39 = DmInv[0]*t31*2.0;
    double t40 = DmInv[2]*t32*2.0;
    double t41 = DmInv[1]*t29*2.0;
    double t42 = DmInv[3]*t30*2.0;
    double t43 = DmInv[1]*t31*2.0;
    double t44 = DmInv[3]*t32*2.0;
    double t45 = t2*t29*2.0;
    double t46 = t3*t30*2.0;
    double t47 = t2*t31*2.0;
    double t48 = t3*t32*2.0;
    double t49 = t37+t38;
    double t50 = t39+t40;
    double t51 = t41+t42;
    double t52 = t43+t44;
    double t57 = t45+t46;
    double t58 = t47+t48;
    double t61 = t33+t34+t35+t36;
    double t53 = t49*t49;
    double t54 = t50*t50;
    double t55 = t51*t51;
    double t56 = t52*t52;
    double t59 = t57*t57;
    double t60 = t58*t58;
    double t62 = 1.0/t61;
    double t63 = sqrt(t61);
    double t64 = 1.0/t63;
    double t65 = t28+t63;
    double t66 = (t6*t53*t62)/4.0;
    double t67 = (t7*t54*t62)/4.0;
    double t68 = (t8*t55*t62)/4.0;
    double t69 = (t9*t56*t62)/4.0;
    double t70 = (t10*t59*t62)/4.0;
    double t71 = (t11*t60*t62)/4.0;
    double t72 = t27+t66+t67+t68+t69+t70+t71;
    double t73 = 1.0/t72;
    deltax[0] = t6*t49*t64*t65*t73*(-1.0/2.0);
    deltax[1] = t7*t50*t64*t65*t73*(-1.0/2.0);
    deltax[2] = t8*t51*t64*t65*t73*(-1.0/2.0);
    deltax[3] = t9*t52*t64*t65*t73*(-1.0/2.0);
    deltax[4] = (t10*t57*t64*t65*t73)/2.0;
    deltax[5] = (t11*t58*t64*t65*t73)/2.0;
    Minv[0] = t6;
    Minv[1] = t7;
    Minv[2] = t8;
    Minv[3] = t9;
    Minv[4] = t10;
    Minv[5] = t11;
    gradC[0] = (t49*t64)/2.0;
    gradC[1] = (t50*t64)/2.0;
    gradC[2] = (t51*t64)/2.0;
    gradC[3] = (t52*t64)/2.0;
    gradC[4] = t57*t64*(-1.0/2.0);
    gradC[5] = t58*t64*(-1.0/2.0);
    deltaLambda[0] = -t65*t73;
    alphaTilde[0] = t27;
    constraintC[0] = t63;
}
