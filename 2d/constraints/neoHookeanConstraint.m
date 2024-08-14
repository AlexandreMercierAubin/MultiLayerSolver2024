classdef neoHookeanConstraint

    methods(Static)
        function [deltaLambda, deltax, forceStackH, alphaTilde, constraintC, gradC, Minv]=project(xi, DmInv, lambda, materialMu, materialLambda, mass, area, h)
            %This is incomplete, don't use it, instead go for the mex
            %neoHookean constraints

            %deformed material frame
            p0 = xi(1:2);
            p1 = xi(3:4);
            p2 = xi(5:6);

            %computing deformation gradient
            Ds  = [p0-p2,p1-p2];
            F = Ds*DmInv;

            %stvk energy
            FtF = F'*F;
            gamma = 1 + (materialMu/materialLambda);
            C_H = det(F)-gamma;
            C_H2 = C_H*C_H;
            alphaH = 1/(area*materialLambda);

            C_D = FtF(1,1)+FtF(2,2)-2;
            alphaD = 1/(area*materialMu);

            psiH = 0.5 *materialLambda *  C_H2;
            psiD = 0.5 * materialMu * C_D;

            constraintH = area * psiH;
            constraintD = area * psiD;

            piolaStress = F*(materialMu*C_H2 + materialLambda*C_D);
            dpsidF = piolaStress * DmInv';
            forceStackH = -area * dpsidF;


            alphaTilde = alpha /(h*h);

            force3 = -forceStackH(:,1) - forceStackH(:,2);
            forceStackH = [force3,forceStackH];
            gradC = -reshape(forceStackH,[], 1);

            numerator = (-constraintC - alphaTilde*lambda);
            Minv = diag(1./mass);
            denumerator = gradC'*Minv*gradC + alphaTilde;
            deltaLambda = numerator/denumerator;
            deltax = Minv*gradC.*deltaLambda;
        end
    end
end

