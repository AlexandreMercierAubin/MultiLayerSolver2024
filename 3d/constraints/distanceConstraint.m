classdef distanceConstraint

    methods(Static)
        function [deltaLambda, deltax1, deltax2]=project(h, p1, p2,mass1,mass2, lambda, alpha, distance)
            assert(all(size(p1) == size(p2)));
            minv1 = 1./mass1;
            minv2 = 1./mass2;
            minvSum = minv1+minv2;

            direction = p1-p2;
            vectorLength = vecnorm(direction,2,2);
            normalized = direction./vectorLength;
            normalized(:,:,vectorLength == 0) = 0; %prevents nan when distance = 0
            constraintC = abs(vectorLength - distance);
            
            h2 = h*h;
            alphaTilde = alpha./h2;

            numerator = -(constraintC + alphaTilde.*lambda);
            denumerator = minvSum + alphaTilde;
            deltaLambda = numerator./denumerator;

            impulse = pagemtimes(normalized,deltaLambda);
            deltax1 = minv1.*impulse;
            deltax2 = -minv2.*impulse;
        end
    end
end

