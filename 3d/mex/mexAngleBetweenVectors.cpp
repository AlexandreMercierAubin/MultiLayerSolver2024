/*==========================================================
 * computes the angle between two vectors wrt a rotation axis
 * N1 is n by 3 normals
 * N2 is n by 3 normals
 * RotationAxis is a N by 3 vector
 * To compile type: 
 mex -R2018a mexAngleBetweenVectors.cpp  -I"../../lib/gptoolbox/mex/external/libigl/external/eigen" -I"../../lib/gptoolbox/mex/external/libigl/include"

 *========================================================*/

#include "mex.h"
#include <math.h>

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    Eigen::MatrixXd N1, N2, RotationAxis;
	
	igl::matlab::parse_rhs_double(prhs+0,N1);
	igl::matlab::parse_rhs_double(prhs+1,N2);
	igl::matlab::parse_rhs_double(prhs+2,RotationAxis);
    const int n = N1.rows();
    Eigen::MatrixXd theta(n,1);
	
    for(size_t i = 0; i < n; ++i){
		const Eigen::Vector3d n1 = N1.row(i);
		const Eigen::Vector3d n2 = N2.row(i);
		const Eigen::Vector3d r = RotationAxis.row(i);
		
		const Eigen::Vector3d ncross = r.cross(n1);
		const Eigen::Vector3d ncrossnorm = ncross/ncross.norm();
		
		Eigen::Matrix3d E;
		E.row(0) = r;
		E.row(1) = n1;
		E.row(2) = ncrossnorm;
		const Eigen::Vector3d c = E * n2;
		theta(i,0) = atan2(c[2],c[1]);
    }
	
	igl::matlab::prepare_lhs_double(theta,plhs+0);
}
