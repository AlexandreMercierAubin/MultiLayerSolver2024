function [deltax, gradC, deltaLambda, constraintC] = matStvkVoigtConstraint3D(x, xprev, DmInv, alphaTilde, lambda, minv, volume, betaTilde, h)%#codegen
    %beta is a rayleigh damping term as a damping stiffness matrix

    [elemGradC,elementConstraintsC] = StVKconstraintGrad3D(DmInv,x);
    % [gradC,constraintC] = StVKconstraintGrad3D_mex(DmInv,x);
    % [gradC,constraintC] = StVKAnalyticalConstraintGrad(DmInv,x); TODO
    constraintC = elementConstraintsC.*volume;
    gradC = elemGradC.*volume;
    
    gradCM = gradC.*pagetranspose(repmat(minv,[1,6,1]));
    
    gamma = pagemtimes(alphaTilde,betaTilde)./h;
    globaldx = x-xprev;
    
    numGamma = pagemtimes(pagemtimes(gamma,gradC),globaldx);
    DenomGamma = (eye(6,6)+gamma);

    % deltaLambda = justSolve(constraintC,alphaTilde,lambda,numGamma, DenomGamma, gradC,gradCM);
    deltaLambda = decomposeSolve(constraintC,alphaTilde,lambda,numGamma, DenomGamma, gradC,gradCM);

    deltax = pagemtimes(gradCM,'transpose',deltaLambda,'none');
end

function deltaLambda = justSolve(constraintC,alphaTilde,lambda,numGamma, DenomGamma, gradC,gradCM)
    numerator = -constraintC - pagemtimes(alphaTilde,lambda) - numGamma;
    denominator = pagemtimes(DenomGamma, pagemtimes(gradCM,'none',gradC,'transpose')) + alphaTilde;

    deltaLambda = pagemldivide(denominator,numerator);
end

function deltaLambda = decomposeSolve(constraintC,alphaTilde,lambda,numGamma, DenomGamma, gradC,gradCM)
    topNumerator = -constraintC(1:3,:,:) - pagemtimes(alphaTilde(1:3,1:3,:),lambda(1:3,:,:)) - numGamma(1:3,:,:);
    topDenominator = pagemtimes(DenomGamma(1:3,1:3,:), pagemtimes(gradCM(1:3,:,:),'none',gradC(1:3,:,:),'transpose')) + alphaTilde(1:3,1:3,:);
    
    numerator = -constraintC(4:6,:,:) - alphaTilde(4,4,:).*lambda(4:6,:,:) - numGamma(4:6,:,:);
    bottom1 = pagemtimes(gradCM(4,:,:),'none',gradC(4,:,:),'transpose');
    bottom2 = pagemtimes(gradCM(5,:,:),'none',gradC(5,:,:),'transpose');
    bottom3 = pagemtimes(gradCM(6,:,:),'none',gradC(6,:,:),'transpose');
    bottom = cat(1,bottom1,bottom2,bottom3);
    denominator = DenomGamma(4,4,:).*bottom + alphaTilde(4,4,:);

    topDeltaLambda = pagemldivide(topDenominator,topNumerator);

    deltaLambda = cat(1,topDeltaLambda, numerator./denominator);
end

