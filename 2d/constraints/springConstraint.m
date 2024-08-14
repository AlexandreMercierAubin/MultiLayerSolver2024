classdef springConstraint

    methods(Static)
        function [constraintC, alphaTilde]=evaluate(h, p1, p2, alpha, restLength, stiffness)
            direction = p2-p1;
            vectorNorm = sqrt(sum(direction.^2));

            springExtension = (vectorNorm - restLength);
            constraintC = 0.5*stiffness*springExtension^2;%hookean spring

            h2 = h*h;
            alphaTilde = alpha/h2;
        end

        function [deltaLambda, deltax1, deltax2,constraintC, alphaTilde]=project(h, p1, p2, w1, w2, mass1, mass2, lambda, alpha, beta, restLength, stiffness)
            minv1 = 1/mass1;
            minv2 = 1/mass2;
            minvSum = w1+w2;

            direction = p2-p1;
            vectorNorm = sqrt(sum(direction.^2));
            if vectorNorm == 0
                normalized = [1;0];
            else
                normalized = direction./vectorNorm;
            end

            springExtension = (vectorNorm - restLength);
            % constraintC = 0.5*stiffness*springExtension^2;%hookean spring
            constraintC = 0.5*stiffness*springExtension; %StVK spring
            if isnan(constraintC)
                error = 1
            end

            h2 = h*h;
            alphaTilde = alpha/h2;
            
            gradConstraint = stiffness*normalized*springExtension;

            numerator = -(constraintC + alphaTilde*lambda);
            denumerator = minvSum + alphaTilde;
            deltaLambda = numerator/denumerator;
            deltax1 = -minv1*deltaLambda*gradConstraint;
            deltax2 = minv2*deltaLambda*gradConstraint;
            if any(isnan(deltax1)) || isnan(deltaLambda)
                error = 1
            end
        end
    end
end

