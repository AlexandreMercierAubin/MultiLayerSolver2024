 classdef rigidBoundaryConstraint2D

    methods(Static)
        function [deltaLambda, deltaxElastic, deltaCOM, deltaOmega]=projectVectorized(h, pElastic, pRigid, relativeRigidPos, minvElastic, lambda, alpha, rigidBodyMass, iterateRigidCOM, iterateRigidRotation, iterateInertia, nb)
            direction = pRigid-pElastic;
            constraintC = vecnorm(direction,2,2);

            n = direction./constraintC;
            n(:,:,constraintC==0)=0; %could happen when the constraint is satisfied

            Iinv = 1./iterateInertia;
            invMassRigid = 1./rigidBodyMass;
            invMassRigidScaled = invMassRigid./nb; %it says split mass and inertia in the mass splitting paper
 
            crossrin = cross2D(relativeRigidPos,n);
            rigidBodyKineticAtPoint = Iinv .* crossrin.^2;
            w = invMassRigidScaled + rigidBodyKineticAtPoint;

            minvSum = minvElastic+w;
            
            h2 = h*h;
            alphaTilde = alpha/h2;

            numerator = -constraintC - alphaTilde.*lambda;
            denominator = minvSum + alphaTilde;
            deltaLambda = numerator./denominator;

            p = reshape(deltaLambda.*n,1,2,[]);
            deltaxElastic = minvElastic.*p;

            deltaCOM = (p.*invMassRigid);
            deltaOmega = (0.5*Iinv).*cross2D(relativeRigidPos,p);
        end

        function [deltaLambda, deltaxElastic, rigidCOM, rigidRotation]=project(h, pElastic, pRigid, relativeRigidPos, massElastic, lambda, alpha, rigidBodyMass, iterateRigidCOM, iterateRigidRotation, iterateInertia, nb)
            assert(size(pRigid,1)==1)
            assert(size(pElastic,1)==1)
            %for the mass splitting term nb, we divide the mass of the
            %constraint by the number of constraints on the rigid body, but
            %let the impulse use the full mass.
            minv = 1/massElastic;

            direction = pElastic-pRigid;
            constraintC = norm(direction);
            if constraintC == 0
                %no need to solve, the constraint is exactly satisfied
                rigidRotation = iterateRigidRotation;
                deltaLambda = 0;
                deltaxElastic = zeros(size(pElastic));
                rigidCOM = iterateRigidCOM;
                return;
            end
            n = direction./constraintC;

            Iinv = 1/iterateInertia;
            invMassRigid = 1/rigidBodyMass;
            invMassRigidScaled = invMassRigid/nb; %it says split mass and inertia in the mass splitting paper
 
            rigidBodyKineticAtPoint = Iinv * cross2D(relativeRigidPos,n)^2;
            w = invMassRigidScaled + rigidBodyKineticAtPoint;

            minvSum = minv+w;
            
            h2 = h*h;
            alphaTilde = alpha/h2;

            numerator = -constraintC - alphaTilde*lambda;
            denominator = minvSum + alphaTilde;
            deltaLambda = numerator/denominator;

            p = deltaLambda*n;
            deltaxElastic = minv*p;
            deltaCOM = (p.*invMassRigid);
            rigidCOM = iterateRigidCOM - deltaCOM;

            deltaOmega = 0.5*Iinv*cross2D(relativeRigidPos,p);
            rigidRotation = iterateRigidRotation - deltaOmega;
        end
    end
end

