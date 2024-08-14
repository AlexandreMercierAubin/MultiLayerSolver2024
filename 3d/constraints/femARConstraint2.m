classdef femARConstraint2

    methods(Static)
        function [deltaLambda, deltax, forceStackH, w, positionalImpulse]=project(xi, mesh3d, alpha, beta, lambda, oldp, h, isElasticDof, rigidBodyIDs, tet, mass, r)
            %fetching element material properties
            materialMu = mesh3d.elMu(tet);
            materialLambda = mesh3d.elLambda(tet);
            
            %deformed material frame
            p0 = xi(1:3)';
            p1 = xi(4:6)';
            p2 = xi(7:9)';
            p3 = xi(10:12)';

            %computing deformation gradient
            Ds  = [p1-p0;p2-p0;p3-p0]';
            DmInv = mesh3d.DmInv{tet};
            F = Ds*DmInv;
%             F = reshape(mesh3d.localB{tet} * xi,3,3);

            %green strain
            E = 0.5 *( F' * F - eye(3) );

            %stvk energy
            froE2 = sum( E.*E, 'all');
            trE = trace(E);
            psi = materialMu * froE2 + materialLambda * 0.5 * trE*trE;

            constraintC = mesh3d.elV(tet) * psi;

            piolaStress = F*(2*materialMu*E + materialLambda*trace(E)*eye(3));
            dpsidF = piolaStress * DmInv';
            forceStackH = -mesh3d.elV(tet) * dpsidF;

            alphaTilde = alpha /(h*h);
            betaTilde = beta*(h*h);
            gamma = alphaTilde *betaTilde / h;

            displacement = xi-oldp;
            force4 = -forceStackH(:,1) - forceStackH(:,2) - forceStackH(:,3);
            forceStackH = [force4,forceStackH];
            gradC = reshape(forceStackH,12, []);
            gradDisplacement = gradC'*displacement;

            minv = zeros(12,1);
            minv(isElasticDof) = 1./mass(isElasticDof);
            isRigidDof = ~isElasticDof;
            rigidStack = gradC(isRigidDof);
            w = zeros(numel(rigidStack),1);
            
            for i = 1:numel(rigidBodyIDs)
                body = mesh3d.RigidBodies(rigidBodyIDs(i));
                n = normalize(-rigidStack(i*3-2:i*3))';
                w(i*3-2:i*3) = rigidGeneralizedInverseMass(1/body.Mass, r(i,:), body.Inertia, n);
            end
%             minv(isRigidDof) = w;
            minv = 1./mass;
            Minv = diag(minv);

            numerator = (-constraintC - alphaTilde*lambda - gamma*gradDisplacement);
            denumerator = (1+gamma)*gradC'*Minv*gradC + alphaTilde;
            deltaLambda = numerator/denumerator;
            deltax = -Minv*gradC.*deltaLambda;

            positionalImpulse = zeros(numel(rigidStack),1);
            for i = 1:numel(rigidBodyIDs)
                n = normalize(-rigidStack(i*3-2:i*3))';
                positionalImpulse(i*3-2:i*3) = n*deltaLambda;
            end
        end
    end
end

