classdef femARConstraint

    methods(Static)
        function [deltaLambda, deltax, forceStackH]=project(cache, xi, mesh3d, alpha, beta, impulseLambda, tet,tetDofsIds,tetDofsIdsVector, G, xdoti, h)
            %fetching element material properties
            materialMu = mesh3d.elMu(tet);
            materialLambda = mesh3d.elLambda(tet);

            %fetching vertices of tet
            dofs = (xi(tetDofsIds));
            
            %deformed material frame
            p0 = dofs(1,:);
            p1 = dofs(2,:);
            p2 = dofs(3,:);
            p3 = dofs(4,:);

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

            xTilde = cache.oldp();
            displacement = xi(tetDofsIdsVector)-cache.oldp(tetDofsIdsVector);
            force4 = -forceStackH(:,1) - forceStackH(:,2) - forceStackH(:,3);
            forceStackH = [force4,forceStackH];
            gradC = reshape(forceStackH,12, []);

            M = diag(mesh3d.mass(tetDofsIdsVector));
            Minv = diag(1./mesh3d.mass(tetDofsIdsVector));
            b = M*(displacement - h*G*xdoti);

            
            smalldeltaxk = xTilde - xi;
            smalldeltaxk =smalldeltaxk(tetDofsIdsVector);%TODO: validate. Might not work
            smalldeltalambdak = Minv*gradC*impulseLambda;
            numerator = -(gamma*gradC'*(smalldeltalambdak + smalldeltaxk) + constraintC + alphaTilde*impulseLambda);
            denumerator = (1+gamma)*gradC'*Minv*gradC + alphaTilde;

            deltaLambda = numerator/denumerator;
            deltax = -Minv*(b + gradC.*deltaLambda);
        end
    end
end
