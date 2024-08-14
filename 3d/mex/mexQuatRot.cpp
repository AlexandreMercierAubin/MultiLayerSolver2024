////////////////////////////////////////////////////////////
// 
// Name:   qrot3d.cpp
//
// Author: Steven Michael
//         (smichael@ll.mit.edu)
// 
//
// Date:   2023/01/12
//
// Description:
//
//   The following is a compiled MATLAB function piece to do 
//   quaternion rotation  in 3 dimensions
//   Read the corresponding .m file for usage.
//   Works with both single & double precision
//
////////////////////////////////////////////////////////////

#include <mex.h>
#include <string.h>
#include <math.h>

#ifndef QUAT_UTILS
#define QUAT_UTILS

typedef double DQuat[4];
typedef double DVec[3];

#ifndef min
#define min(a,b) (a < b ? a : b)
#endif


static void dqmult(DQuat a, DQuat b, DQuat c)
{
  c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  c[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
  c[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
  c[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
  return;
} // end of dqmult

#endif

// A function that rotates a single vector by the desired
// quaternion
static void rotate_double(DVec vout,DVec vin,DQuat q, int n)
{
  // Construct the quaternion conjugate
  DQuat qc;
  int i;
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];

  // Do the quaternion multiplitcation
  // for all the input vectors
  for(i=0;i<n;i++) {
    DQuat qv,qt,qf;
    qv[0] = 0;
    qv[1] = vin[i];
    qv[2] = vin[i+n];
    qv[3] = vin[i+2*n];
    
    dqmult(q,qv,qt);
    dqmult(qt,qc,qf);
    vout[i] = qf[1];
    vout[i+n] = qf[2];
    vout[i+2*n] = qf[3];
  }
  return;
} // end of rotate_double

// Construct a quaternion from a rotation vector
// and angle

void d_construct_quaternion(DVec w, double theta, DQuat q)
{
  // Normalize w
  double norm = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  double ct2,st2;
  if(norm == 0) {
    mexErrMsgTxt("Norm of w cannot be zero\n");
    return;
  }
  w[0] = w[0]/norm;
  w[1] = w[1]/norm;
  w[2] = w[2]/norm;

  // Construct quaternions
  ct2 = cos(theta/2.0f);
  st2 = sin(theta/2.0f);
  q[0] = ct2;
  q[1] = st2*w[0];
  q[2] = st2*w[1];
  q[3] = st2*w[2];  
} // end of d_construct_quaternion

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // The rotation angle, in double
    // or single float format
    if(mxGetNumberOfElements(prhs[1])!=4) {
      mexErrMsgTxt("2nd input must have length 4");  
    }
    DQuat  dq;
    int    ndims;
    mwSize dims[10];
    
    double *d = (double *)mxGetPr(prhs[1]);
    dq[0] = d[0];dq[1] = d[1];dq[2] = d[2];dq[3] = d[3];
    
    
    // Get the dimensions of the input
    ndims = mxGetNumberOfDimensions(prhs[0]);
    
    memcpy(dims,mxGetDimensions(prhs[0]),
     sizeof(mwSize)*ndims);

    // Make sure the dimensions are what they
    // are supposed to be
    if(dims[1] != 3) {
    mexErrMsgTxt("Input Vectors Must Be (N x 3)\n");
    return;
    }

    double *din;
    double *dout;
    plhs[0] = 
    mxCreateNumericMatrix(dims[0],dims[1],
            mxDOUBLE_CLASS,
            mxREAL);
    din = (double *)mxGetPr(prhs[0]);
    dout = (double *)mxGetPr(plhs[0]);
    rotate_double(dout,din,dq,dims[0]);

    return;
}

