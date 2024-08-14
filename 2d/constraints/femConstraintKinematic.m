classdef femConstraintKinematic

    methods(Static)
        function [alphaTilde, constraintC, gradC]=project(cache, positions, mesh2d, alpha, tri,h)
            %fetching element material properties
            materialMu = mesh2d.elMu(tri);
            materialLambda = mesh2d.elLambda(tri);
            
            %deformed material frame
            p0 = positions(1,:);
            p1 = positions(2,:);
            p2 = positions(3,:);

            %computing deformation gradient
            Ds  = [p0-p2;p1-p2]';

            DmInv = mesh2d.DmInv(:,:,tri);
            F = Ds*DmInv;

            %green strain
            E = 0.5 *( F' * F - eye(2) );

            %stvk energy
            EtE = E'*E;
            froE2 = EtE(1,1) + EtE(2,2);
            trE = E(1,1)+E(2,2);
            psi = materialMu * froE2 + materialLambda * 0.5 * trE*trE;

            constraintC = mesh2d.elA(tri) * psi;

            %page 18 femdefo. This is dpsidF
            piolaStress = F*(2*materialMu*E + materialLambda*trE*eye(2)); 
            %page 29 forces = -area*dpsidx
            dpsidx = piolaStress * DmInv';
            W = -mesh2d.elA(tri) * dpsidx;
            force3 = -W(:,1) - W(:,2);
            forceStackH = [W,force3];

            %so logically gradC = -forces because f = -area grad psi
            gradC = -reshape(forceStackH,[], 1);

            h2 = h*h;
            alphaTilde = alpha /h2;
        end
    end
end

