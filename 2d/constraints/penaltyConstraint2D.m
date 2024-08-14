classdef penaltyConstraint2D

    methods(Static)
        function [deltax, deltaLambda, C]=projectElastic(x, lambda, minv, contactPoint, contactNormal, contactResistance, alpha, h)
            %gradC is the contact normal
            contactVector = x - contactPoint;
            phi = dot(contactVector, contactNormal,2);
            phi(phi>0) = 0;% no interpenetration

            C = contactResistance .*phi;
            if numel(contactPoint) == 0
                deltax = zeros(size(x));
                deltaLambda = zeros(size(lambda));
                return;
            end

            h2 = h*h;
            alphaTilde = alpha/h2;

            numerator = -(C + alphaTilde.*lambda);
            gradCminv = contactNormal.*2*minv;
            denumerator = pagemtimes(gradCminv,'none',contactNormal,'transpose') + alphaTilde;
            deltaLambda = numerator./denumerator;

            deltax = pagetranspose(pagemtimes(gradCminv,'transpose',deltaLambda,'none'));
        end

        function [deltaLambda, deltaCOM, deltaOmega, C]=projectRigid(h, x, contactPoint, contactNormal, contactResistance, relativeRigidPos, lambda, alpha,  rigidBodyMass, minv, iterateInertia, nb)

            contactVector = x - contactPoint;
            phi = dot(contactVector, contactNormal,2);
            phi(phi>0) = 0;% no interpenetration
            
            C = contactResistance .* phi;
            if numel(contactPoint) == 0
                deltaLambda = [];
                deltaCOM = [];
                deltaOmega = [];
                return;
            end

            Iinv = 1./iterateInertia;
            invMassRigid = 1./rigidBodyMass;
            invMassRigidScaled = invMassRigid./nb;

            crossrin = cross2D(relativeRigidPos,contactNormal);
            rigidBodyKineticAtPoint = Iinv .* crossrin.^2;
            w = invMassRigidScaled + rigidBodyKineticAtPoint;
            
            h2 = h*h;
            alphaTilde = alpha./h2;

            numerator = -C - alphaTilde.*lambda;
            denumerator = 2 * w + alphaTilde;
            deltaLambda = numerator./denumerator;

            p = deltaLambda.*contactNormal;
            deltaCOM = (p.*invMassRigid);
            deltaOmega = (0.5*Iinv).*cross2D(relativeRigidPos,p);
        end
    end
end

