classdef rigidBoundaryConstraint

    methods(Static)
        function [deltaLambda, deltaxElastic, deltaCOM, deltaOmega]=project(h, pElastic, pRigid, relativeRigidPos, massInvElastic, lambda, alphaTilde, rigidBodyMass, iterateInertia, nb)
            assert(size(pElastic,3) == size(pRigid,3));

            direction = pRigid-pElastic;
            constraintC = vecnorm(direction,2,2);

            if numel(rigidBodyMass) == 0
                deltaLambda = [];
                deltaCOM = [];
                deltaOmega = [];
                deltaxElastic = zeros(size(pElastic));
                return;
            end
            n = direction./constraintC;
            n(isnan(n(:))) = 0; %makes sure we don't have issues when the constraint is perfectly satisfied

            invMassRigid = 1./rigidBodyMass;
            % invMassRigidScaled = invMassRigid./nb;
            w = rigidGeneralizedInverseMass(invMassRigid, relativeRigidPos, iterateInertia, n);

            minvSum = massInvElastic+w;

            numerator = -constraintC - alphaTilde.*lambda;
            denumerator = minvSum + alphaTilde;
            deltaLambda = numerator./denumerator;

            p = deltaLambda.*n;
            deltaCOM = p.*invMassRigid;
            deltaxElastic = -massInvElastic.*p;

            deltaOmega = cross(relativeRigidPos,p./h,2);
        end
    end
end

