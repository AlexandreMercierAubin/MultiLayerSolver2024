function R = expRodrigues( omega, h )
    %EXPRODRIGUES returns a rotation matrix exp(h * omegahat);
    % omega     3D angular velocity
    % h         time step
    assert(size(omega,1) == 1);
    
    mag = vecnorm( omega, 2);
    w = omega ./ mag;
    w(isnan(w))=0; %can happen when there is no rotation
    
    theta = h .* mag;

    K = crossProductMatrix(w);
    R = eye(3) + sin(theta).*K + ((1-cos(theta))).*pagemtimes(K,K);
end

%Legacy code for good measure
% c = cos(theta);
% s = sin(theta);
% 
% c1 = 1 - c;
% R = zeros(3,3);
% R(1,1) =  c        + w(1) * w(1) * c1;
% R(2,1) =  w(3) * s + w(1) * w(2) * c1;
% R(3,1) = -w(2) * s + w(1) * w(3) * c1;
% 
% R(1,2) = -w(3) * s + w(1) * w(2) * c1;
% R(2,2) =  c        + w(2) * w(2) * c1;
% R(3,2) =  w(1) * s + w(2) * w(3) * c1;
% 
% R(1,3) =  w(2) * s + w(1) * w(3) * c1;
% R(2,3) = -w(1) * s + w(2) * w(3) * c1;
% R(3,3) =  c        + w(3) * w(3) * c1;