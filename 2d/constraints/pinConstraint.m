classdef pinConstraint

    methods(Static)
        function [deltalambda, deltax, constraintC]=project(p, pinPos,mass)
            l = pinPos-p;
            constraintC = norm(l);
            deltax = l;
            deltalambda = -constraintC*mass;
        end
    end
end

