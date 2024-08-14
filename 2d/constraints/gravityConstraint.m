classdef gravityConstraint

    methods(Static)
        function [deltax, deltaLambda, C]=project(x, xprev, lambda, mass, alpha, g, h)
            dx = x - xprev;
            gravityDx = g*h*h;
            gradC= -max(dx - gravityDx ,0);
            C = -gradC;

            h2 = h*h;
            alphaTilde = alpha/h2;

            numerator = -(C + alphaTilde*lambda);
            minv = (1./mass);
            denominator = (gradC.*minv).*gradC + alphaTilde;
            deltaLambda = numerator./denominator;

            deltax = -gradC.*deltaLambda;
        end
    end
end

