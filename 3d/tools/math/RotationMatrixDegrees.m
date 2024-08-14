function [R] = RotationMatrixDegrees(degrees)
% RotationMatrixDegrees Computes a rotation matrix from Euler angles
%
%  [R] = RotationMatrixDegrees(degrees)

    R = eul2rotm(toRadians("degrees",degrees),'XYZ');
end

