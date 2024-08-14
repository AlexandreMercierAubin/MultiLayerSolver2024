% mexEdiffNorm3D.cpp - computes 3D strain difference Frobeneneous norm squared
%
% [ eDotNormsSquared ] = mexEdotNorm2D( Fa, Fb );
%
% Input:
%  Fa           9x#E deformation gradient of each element (or a vector of size 9x#E)
%  Fb           9x#E deformation gradient of each element (or a vector of size 9x#E)
%  h            timestep
%
% Output:
%  EdiffNormSquared        #E Frobeneus norm squared of (strain difference
%                             divided by h), i.e., an approximation of strain rate
%
% To compile type: mex -R2018a mexEdiffNorm3D.cpp