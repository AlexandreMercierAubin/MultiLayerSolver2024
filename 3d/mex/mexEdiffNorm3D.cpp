/*==========================================================
 * mexEdiffNorm.cpp - computes 3D strain difference Frobeneneous norm squared
 *
 * [ eDotNormsSquared ] = mexEdotNorm( Fa, Fb, h );
 *
 * Input:
 *  Fa      9x#E deformation gradient of each element (or a vector of size 9x#E) at one time step
 *  Fb      9x#E deformation gradient of each element (or a vector of size 9x#E) at another time step
 *  h       timestep
 *
 * Output:
 *  EdiffNormSquared     #E Frobeneus norm squared of strain difference 
 *                          divided by h), i.e., an approximation of strain rate
 *
 * To compile type: mex -R2018a mexEdiffNorm3D.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h>

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *Fa = mxGetDoubles(prhs[0]);
    size_t sizeF = mxGetM(prhs[0]);
    double *Fb = mxGetDoubles(prhs[1]);
    double h = mxGetScalar( prhs[2] );

    int sizeEdiffNorm = sizeF/9.0;
    plhs[0] = mxCreateDoubleMatrix( sizeEdiffNorm, 1, mxREAL );
    double *EdiffNormSquared = mxGetDoubles(plhs[0]);
    
    for(size_t i = 0; i < sizeEdiffNorm; ++i) {
        int pos = i*9;
        double Fa1_1 = Fa[pos];
        double Fa2_1 = Fa[pos+1];
        double Fa3_1 = Fa[pos+2];
        double Fa1_2 = Fa[pos+3];
        double Fa2_2 = Fa[pos+4];
        double Fa3_2 = Fa[pos+5];
        double Fa1_3 = Fa[pos+6];
        double Fa2_3 = Fa[pos+7];
        double Fa3_3 = Fa[pos+8];
        
        double Fb1_1 = Fb[pos];
        double Fb2_1 = Fb[pos+1];
        double Fb3_1 = Fb[pos+2];
        double Fb1_2 = Fb[pos+3];
        double Fb2_2 = Fb[pos+4]; 
        double Fb3_2 = Fb[pos+5]; 
        double Fb1_3 = Fb[pos+6]; 
        double Fb2_3 = Fb[pos+7]; 
        double Fb3_3 = Fb[pos+8]; 
        
        double t2 = 1.0/(h*h);
        double t0 = t2*pow((Fa1_1*Fa1_2)/2.0+(Fa2_1*Fa2_2)/2.0+(Fa3_1*Fa3_2)/2.0-(Fb1_1*Fb1_2)/2.0-(Fb2_1*Fb2_2)/2.0-(Fb3_1*Fb3_2)/2.0,2.0)*2.0+t2*pow((Fa1_1*Fa1_3)/2.0+(Fa2_1*Fa2_3)/2.0+(Fa3_1*Fa3_3)/2.0-(Fb1_1*Fb1_3)/2.0-(Fb2_1*Fb2_3)/2.0-(Fb3_1*Fb3_3)/2.0,2.0)*2.0+t2*pow((Fa1_2*Fa1_3)/2.0+(Fa2_2*Fa2_3)/2.0+(Fa3_2*Fa3_3)/2.0-(Fb1_2*Fb1_3)/2.0-(Fb2_2*Fb2_3)/2.0-(Fb3_2*Fb3_3)/2.0,2.0)*2.0+t2*pow((Fa1_1*Fa1_1)/2.0+(Fa2_1*Fa2_1)/2.0+(Fa3_1*Fa3_1)/2.0-(Fb1_1*Fb1_1)/2.0-(Fb2_1*Fb2_1)/2.0-(Fb3_1*Fb3_1)/2.0,2.0)+t2*pow((Fa1_2*Fa1_2)/2.0+(Fa2_2*Fa2_2)/2.0+(Fa3_2*Fa3_2)/2.0-(Fb1_2*Fb1_2)/2.0-(Fb2_2*Fb2_2)/2.0-(Fb3_2*Fb3_2)/2.0,2.0)+t2*pow((Fa1_3*Fa1_3)/2.0+(Fa2_3*Fa2_3)/2.0+(Fa3_3*Fa3_3)/2.0-(Fb1_3*Fb1_3)/2.0-(Fb2_3*Fb2_3)/2.0-(Fb3_3*Fb3_3)/2.0,2.0);

        EdiffNormSquared[i] = t0;
    }
}
