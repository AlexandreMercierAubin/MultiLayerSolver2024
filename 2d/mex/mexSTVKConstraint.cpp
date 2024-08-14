/*==========================================================
 * mexSTVKConstraint.cpp 
 *
 * To compile type: mex -R2018a mexSTVKConstraint.cpp
 *========================================================*/

#include "mex.h"
#include <math.h>
#include <cmath>

double cross2D(double v11, double v12, double v21, double v22){
    return v11*v22 - v12*v21;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
        /* check for proper number of arguments */
    if ( nrhs != 11 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nrhs","11 inputs required.");
    }
    if ( nlhs != 6 ) {
        mexErrMsgIdAndTxt("ARP:constraint:nlhs","6 outputs required.");
    }

    // could read like this for cleaner file, but likely a bit slower: double volume = *mxGetDoubles(prhs[6]);

    double *h = mxGetDoubles(prhs[0]);
    double* x = mxGetDoubles(prhs[1]);
    double* xOld = mxGetDoubles(prhs[2]);
    double* DmInv = mxGetDoubles(prhs[3]);
    double *materialMu = mxGetDoubles(prhs[4]);
    double *materialLambda = mxGetDoubles(prhs[5]);
    double *volume = mxGetDoubles(prhs[6]);
    double *alpha = mxGetDoubles(prhs[7]);
    double *beta = mxGetDoubles(prhs[8]);
    double *lambda = mxGetDoubles(prhs[9]);
    double* mass = mxGetDoubles(prhs[10]);

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

    double t2 = volume[0]*volume[0];
    double t3 = DmInv[0]+DmInv[1];
    double t4 = DmInv[2]+DmInv[3];
    double t5 = 1.0/(h[0]*h[0]);
    double t6 = 1.0/(h[0]*h[0]*h[0]*h[0]*h[0]);
    double t7 = 1.0/mass[0];
    double t8 = 1.0/mass[1];
    double t9 = 1.0/mass[2];
    double t10 = 1.0/mass[3];
    double t11 = 1.0/mass[4];
    double t12 = 1.0/mass[5];
    double t13 = -x[4];
    double t14 = -x[5];
    double t15 = -xOld[0];
    double t16 = -xOld[1];
    double t17 = -xOld[2];
    double t18 = -xOld[3];
    double t19 = -xOld[4];
    double t20 = -xOld[5];
    double t21 = alpha[0]*t5;
    double t22 = t13+x[0];
    double t23 = t14+x[1];
    double t24 = t13+x[2];
    double t25 = t14+x[3];
    double t26 = t15+x[0];
    double t27 = t16+x[1];
    double t28 = t17+x[2];
    double t29 = t18+x[3];
    double t30 = t19+x[4];
    double t31 = t20+x[5];
    double t32 = alpha[0]*beta[0]*t6;
    double t33 = lambda[0]*t21;
    double t34 = DmInv[0]*t22;
    double t35 = DmInv[2]*t22;
    double t36 = DmInv[0]*t23;
    double t37 = DmInv[2]*t23;
    double t38 = DmInv[1]*t24;
    double t39 = DmInv[3]*t24;
    double t40 = DmInv[1]*t25;
    double t41 = DmInv[3]*t25;
    double t42 = t32+1.0;
    double t43 = t34+t38;
    double t44 = t35+t39;
    double t45 = t36+t40;
    double t46 = t37+t41;
    double t47 = t43*t43;
    double t48 = t44*t44;
    double t49 = t45*t45;
    double t50 = t46*t46;
    double t51 = DmInv[0]*t43;
    double t52 = DmInv[2]*t44;
    double t53 = DmInv[0]*t45;
    double t54 = DmInv[2]*t46;
    double t55 = DmInv[1]*t43;
    double t56 = DmInv[3]*t44;
    double t57 = DmInv[1]*t45;
    double t58 = DmInv[3]*t46;
    double t59 = t3*t43;
    double t60 = t4*t44;
    double t61 = t3*t45;
    double t62 = t4*t46;
    double t63 = (DmInv[2]*t43)/2.0;
    double t64 = (DmInv[0]*t44)/2.0;
    double t65 = (DmInv[2]*t45)/2.0;
    double t66 = (DmInv[0]*t46)/2.0;
    double t67 = (DmInv[3]*t43)/2.0;
    double t68 = (DmInv[1]*t44)/2.0;
    double t69 = (DmInv[3]*t45)/2.0;
    double t70 = (DmInv[1]*t46)/2.0;
    double t75 = (t4*t43)/2.0;
    double t76 = (t3*t44)/2.0;
    double t77 = (t4*t45)/2.0;
    double t78 = (t3*t46)/2.0;
    double t79 = (t43*t44)/2.0;
    double t80 = (t45*t46)/2.0;
    double t71 = t47/2.0;
    double t72 = t48/2.0;
    double t73 = t49/2.0;
    double t74 = t50/2.0;
    double t81 = t51+t52;
    double t82 = t53+t54;
    double t83 = t55+t56;
    double t84 = t57+t58;
    double t85 = t59+t60;
    double t86 = t61+t62;
    double t87 = t63+t64;
    double t88 = t65+t66;
    double t89 = t67+t68;
    double t90 = t69+t70;
    double t91 = t75+t76;
    double t92 = t77+t78;
    double t109 = t79+t80;
    double t93 = t71+t73-1.0/2.0;
    double t94 = t72+t74-1.0/2.0;
    double t110 = t109*t109;
    double t112 = t71+t72+t73+t74-1.0;
    double t115 = t87*t109*4.0;
    double t116 = t88*t109*4.0;
    double t117 = t89*t109*4.0;
    double t118 = t90*t109*4.0;
    double t119 = t91*t109*4.0;
    double t120 = t92*t109*4.0;
    double t95 = t93*t93;
    double t96 = t94*t94;
    double t97 = t51*t93*2.0;
    double t98 = t53*t93*2.0;
    double t99 = t52*t94*2.0;
    double t100 = t55*t93*2.0;
    double t101 = t54*t94*2.0;
    double t102 = t57*t93*2.0;
    double t103 = t56*t94*2.0;
    double t104 = t58*t94*2.0;
    double t105 = t59*t93*2.0;
    double t106 = t61*t93*2.0;
    double t107 = t60*t94*2.0;
    double t108 = t62*t94*2.0;
    double t111 = t110*2.0;
    double t113 = t112*t112;
    double t121 = materialLambda[0]*t81*t112;
    double t122 = materialLambda[0]*t82*t112;
    double t123 = materialLambda[0]*t83*t112;
    double t124 = materialLambda[0]*t84*t112;
    double t125 = materialLambda[0]*t85*t112;
    double t126 = materialLambda[0]*t86*t112;
    double t114 = (materialLambda[0]*t113)/2.0;
    double t127 = t95+t96+t111;
    double t129 = t97+t99+t115;
    double t130 = t98+t101+t116;
    double t131 = t100+t103+t117;
    double t132 = t102+t104+t118;
    double t137 = t105+t107+t119;
    double t138 = t106+t108+t120;
    double t128 = materialMu[0]*t127;
    double t133 = materialMu[0]*t129;
    double t134 = materialMu[0]*t130;
    double t135 = materialMu[0]*t131;
    double t136 = materialMu[0]*t132;
    double t139 = materialMu[0]*t137;
    double t140 = materialMu[0]*t138;
    double t141 = t114+t128;
    double t143 = t121+t133;
    double t144 = t122+t134;
    double t145 = t123+t135;
    double t146 = t124+t136;
    double t155 = t125+t139;
    double t156 = t126+t140;
    double t142 = t141*volume[0];
    double t147 = t143*t143;
    double t148 = t144*t144;
    double t149 = t145*t145;
    double t150 = t146*t146;
    double t151 = t26*t143*volume[0];
    double t152 = t27*t144*volume[0];
    double t153 = t28*t145*volume[0];
    double t154 = t29*t146*volume[0];
    double t157 = t155*t155;
    double t158 = t156*t156;
    double t163 = t30*t155*volume[0];
    double t164 = t31*t156*volume[0];
    double t159 = t2*t7*t42*t147;
    double t160 = t2*t8*t42*t148;
    double t161 = t2*t9*t42*t149;
    double t162 = t2*t10*t42*t150;
    double t165 = -t163;
    double t166 = -t164;
    double t167 = t2*t11*t42*t157;
    double t168 = t2*t12*t42*t158;
    double t169 = t151+t152+t153+t154+t165+t166;
    double t171 = t21+t159+t160+t161+t162+t167+t168;
    double t170 = t32*t169;
    double t172 = 1.0/t171;
    double t173 = t33+t142+t170;
    deltax[0] = -t7*t143*t172*t173*volume[0];
    deltax[1] = -t8*t144*t172*t173*volume[0];
    deltax[2] = -t9*t145*t172*t173*volume[0];
    deltax[3] = -t10*t146*t172*t173*volume[0];
    deltax[4] = t11*t155*t172*t173*volume[0];
    deltax[5] = t12*t156*t172*t173*volume[0];
    Minv[0] = t7;
    Minv[1] = t8;
    Minv[2] = t9;
    Minv[3] = t10;
    Minv[4] = t11;
    Minv[5] = t12;
    gradC[0] = t143*volume[0];
    gradC[1] = t144*volume[0];
    gradC[2] = t145*volume[0];
    gradC[3] = t146*volume[0];
    gradC[4] = -t155*volume[0];
    gradC[5] = -t156*volume[0];
    deltaLambda[0] = -t172*t173;
    alphaTilde[0] = t21;
    constraintC[0] = t142;
}
