classdef stvkConstraint

    methods(Static)
        function [deltaLambda, deltax, forceStackH, alphaTilde, constraintC, gradC, Minv]=project(oldp, xi, DmInv, alpha, beta, lambda, materialMu, materialLambda, mass, area, h)
            %might want to set alpha with respect to the formulas from A Constraint-based Formulation of Stable Neo-Hookean Materials
            %even if it might be off
            %deformed material frame
            p0 = xi(1:2);
            p1 = xi(3:4);
            p2 = xi(5:6);

            %computing deformation gradient
            Ds  = [p0-p2,p1-p2];
            F = Ds*DmInv;

            %green strain
            E = 0.5 *( F' * F - eye(2) );

            %stvk energy
            EtE = E'*E;
            froE2 = EtE(1,1) + EtE(2,2);
            trE = E(1,1)+E(2,2);
            psi = materialMu * froE2 + materialLambda * 0.5 * trE*trE;

            constraintC = area * psi;
            % constraintC = E(1,1)+E(2,2)+E(1,2);

            %page 18 femdefo. This is dpsidF
            piolaStress = F*(2*materialMu*E + materialLambda*trE*eye(2)); 
            %page 29 forces = -area*dpsidx
            dpsidx = piolaStress * DmInv';
            W = -area * dpsidx;
            force3 = -W(:,1) - W(:,2);
            forceStackH = [W,force3];

            %so logically gradC = -forces because f = -area grad psi
            gradC = -reshape(forceStackH,[], 1);

            h2 = h*h;
            alphaTilde = alpha /h2;
            betaTilde = beta* h2;
            gamma = (alphaTilde *betaTilde) / h;

            displacement = xi-oldp;
            
            
            gradDisplacement = gradC'*displacement;

            numerator = -constraintC - alphaTilde*lambda - gamma*gradDisplacement;
            Minv = diag(1./mass);
            denumerator = (1+gamma)*gradC'*Minv*gradC + alphaTilde;
            deltaLambda = numerator/denumerator;
            deltax = Minv*gradC*deltaLambda;
        end

        function [alphaTilde, constraintC]=evaluate(xi, DmInv, alpha, materialMu, materialLambda, area, h)
            %might want to set alpha with respect to the formulas from A Constraint-based Formulation of Stable Neo-Hookean Materials
            %even if it might be off
            %deformed material frame
            p0 = xi(1:2);
            p1 = xi(3:4);
            p2 = xi(5:6);

            %computing deformation gradient
            Ds  = [p0-p2,p1-p2];
            F = Ds*DmInv;

            %green strain
            E = 0.5 *( F' * F - eye(2) );

            %stvk energy
            EtE = E'*E;
            froE2 = EtE(1,1) + EtE(2,2);
            trE = E(1,1)+E(2,2);
            psi = materialMu * froE2 + materialLambda * 0.5 * trE*trE;

            constraintC = area * psi;

            h2 = h*h;
            alphaTilde = alpha /h2;
        end
    end
end

