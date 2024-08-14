classdef penaltyConstraint

    methods(Static)
        function [deltax, deltaLambda, C]=projectElastic(x, lambda, minv, contactPoint, contactNormal, contactResistance, alpha, h)
            %lambda will be 0 as we are making this a hard constraint
            assert(size(contactPoint,2)==3);

            direction =  x - contactPoint;
            phi = dot(contactNormal, direction,2);
            phi(phi>0) = 0;% no interpenetration
 
            C = contactResistance * phi;
            if numel(contactPoint) == 0
                deltax = zeros(size(x));
                deltaLambda = zeros(size(lambda));
                return;
            end

            h2 = h*h;
            alphaTilde = alpha./h2;

            numerator = -C - pagemtimes(alphaTilde,lambda);
            denumerator = 2*minv(1,1,:) + alphaTilde;

            deltaLambda = numerator./denumerator;
            p = deltaLambda.*contactNormal;
            deltax = minv.*p;
        end

        function [deltaLambda, deltaCOM, deltaOmega]=projectRigid(h, x, contactPoint, contactNormal, contactResistance, relativeRigidPos, lambda, alpha, rigidBodyMass, minv, iterateInertia, nb)

            direction = x - contactPoint;
            phi = dot(contactNormal, direction,2);
            phi(phi>0) = 0; % if positive, then we are not in contact

            if numel(contactPoint) == 0
                deltaLambda = [];
                deltaCOM = [];
                deltaOmega = [];
                return;
            end

            C = contactResistance .* phi;

            invMassRigid = 1./rigidBodyMass;
            %invMassRigidScaled = invMassRigid./nb;
            w = rigidGeneralizedInverseMass(invMassRigid, relativeRigidPos, iterateInertia, contactNormal);
            
            h2 = h*h;
            alphaTilde = alpha./h2;

            numerator = -C - pagemtimes(alphaTilde,lambda);
            %2*w to assume the contact pair has the same mass
            denumerator = 2*w + alphaTilde;
            deltaLambda = numerator./denumerator;

            p = deltaLambda.*contactNormal;
            deltaCOM = p.*invMassRigid;
            
            %the inertia part is added in the xpbdLayer constraint code to
            %reduce the number of calls
            deltaOmega = (1/h).*cross(relativeRigidPos,p,2);
        end
    end
end

