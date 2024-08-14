function [deltax, gradC, deltaLambda, constraintC] = matStvkVoigtConstraint(x, xprev, DmInv, alphaTilde, lambda, minv, betaTilde, h, volume)
    %beta is a rayleigh damping term as a damping stiffness matrix

    % [gradCelement,constraintCelement] = StVKconstraintGrad(DmInv,x);
    [gradCelement,constraintCelement] = StVKAnalyticalConstraintGrad(DmInv,x);
    gradC = volume.*gradCelement;
    constraintC = volume.*constraintCelement;
    
    gradCM = gradC.*pagetranspose(repmat(minv,[1,3,1]));
    

    gamma = pagemtimes(alphaTilde,betaTilde)./h;
    globaldx = x-xprev;
    
    numGamma = pagemtimes(pagemtimes(gamma,gradC),globaldx);
    DenomGamma = (eye(3,3)+gamma);

    % deltaLambda = decomposedSolve(constraintC,alphaTilde,lambda,numGamma, DenomGamma, gradC,gradCM);
    deltaLambda = justSolve(constraintC,alphaTilde,lambda,numGamma, DenomGamma, gradC,gradCM);

    deltax = pagemtimes(gradCM,'transpose',deltaLambda,'none');
end

function deltaLambda = justSolve(constraintC,alphaTilde,lambda,numGamma, DenomGamma, gradC,gradCM)
    numerator = -constraintC - pagemtimes(alphaTilde,lambda) - numGamma;
    denominator = pagemtimes(DenomGamma, pagemtimes(gradCM,'none',gradC,'transpose')) + alphaTilde;

    deltaLambda = pagemldivide(denominator,numerator);
    % topDeltaLambda = zeros(2,1,size(alphaTilde,3));
    % coder.gpu.kernelfun();
    % for page = 1:size(topDeltaLambda,3)
    %     deltaLambda(:,:,page) = denominator(:,:,page) \ numerator(:,:,page);
    % end
end

function deltaLambda = decomposedSolve(constraintC,alphaTilde,lambda,numGamma, DenomGamma, gradC,gradCM)
    topNumerator = -constraintC(1:2,:,:) - pagemtimes(alphaTilde(1:2,1:2,:),lambda(1:2,:,:)) - numGamma(1:2,:,:);
    topDenominator = pagemtimes(DenomGamma(1:2,1:2,:), pagemtimes(gradCM(1:2,:,:),'none',gradC(1:2,:,:),'transpose')) + alphaTilde(1:2,1:2,:);
    
    numerator = -constraintC(3,:,:) - pagemtimes(alphaTilde(3,3,:),lambda(3,:,:)) - numGamma(3,:,:);
    denominator = pagemtimes(DenomGamma(3,3,:),pagemtimes(gradCM(3,:,:),'none',gradC(3,:,:),'transpose')) + alphaTilde(3,3,:);

    topDeltaLambda = pagemldivide(topDenominator,topNumerator);
    % topDeltaLambda = zeros(2,1,size(alphaTilde,3));
    % coder.gpu.kernelfun();
    % for page = 1:size(topDeltaLambda,3)
    %     topDeltaLambda(:,:,page) = topDenominator(:,:,page) \ topNumerator(:,:,page);
    % end

    deltaLambda = cat(1,topDeltaLambda, numerator./denominator);
end

function [gradC,constraintC] = StVKAnalyticalConstraintGrad(DmInv,x)
    p3 = x(5:6,:,:);
    Ds = cat(2,x(1:2,:,:)-p3,x(3:4,:,:)-p3);
    F = pagemtimes(Ds,DmInv);
    FtF = pagemtimes(F,'transpose',F,'none');
    E = 0.5.*(FtF-eye(2));
    constraintC = [E(1,1,:);E(2,2,:);E(1,2,:)];

    bigF = repmat(F,2,1);
    bigDm = DmInv([1,1,2,2],:,:);

    b = pagetranspose(bigDm.*bigF);

    b3 = -(b(:,1:2,:)+b(:,3:4,:));

    permDmInv = DmInv(:,[2,1],:);

    bigF2 = repmat(F,3,1);
    sumDm = sum(permDmInv,1);
    bigPermDm = [permDmInv([1,1],:,:);permDmInv([2,2],:,:);-sumDm([1,1],:,:)];

    gradCShear = pagetranspose(0.5.*sum(bigPermDm.*bigF2,2));

    upperGradC = cat(2,b,b3);
    gradC = cat(1,upperGradC,gradCShear);
end

function [gradC,constraintC] = StVKconstraintGrad(DmInv,x)
    t2 = DmInv(1,1,:)+DmInv(2,1,:);
    t3 = DmInv(1,2,:)+DmInv(2,2,:);
    t4 = -x(5,:,:);
    t5 = -x(6,:,:);
    t6 = t4+x(1,:,:);
    t7 = t5+x(2,:,:);
    t8 = t4+x(3,:,:);
    t9 = t5+x(4,:,:);
    t10 = DmInv(1,1,:).*t6;
    t11 = DmInv(1,2,:).*t6;
    t12 = DmInv(1,1,:).*t7;
    t13 = DmInv(1,2,:).*t7;
    t14 = DmInv(2,1,:).*t8;
    t15 = DmInv(2,2,:).*t8;
    t16 = DmInv(2,1,:).*t9;
    t17 = DmInv(2,2,:).*t9;
    t18 = t10+t14;
    t19 = t11+t15;
    t20 = t12+t16;
    t21 = t13+t17;

    out = [DmInv(1,1,:).*t18;%gradient
            DmInv(1,1,:).*t20;
            DmInv(2,1,:).*t18;
            DmInv(2,1,:).*t20;
            -t2.*t18;
            -t2.*t20;
            DmInv(1,2,:).*t19;
            DmInv(1,2,:).*t21;
            DmInv(2,2,:).*t19;
            DmInv(2,2,:).*t21;
            -t3.*t19;
            -t3.*t21;
            (DmInv(1,1,:).*t19)./2.0+(DmInv(1,2,:).*t18)./2.0;
            (DmInv(1,1,:).*t21)./2.0+(DmInv(1,2,:).*t20)./2.0;
            (DmInv(2,1,:).*t19)./2.0+(DmInv(2,2,:).*t18)./2.0;
            (DmInv(2,1,:).*t21)./2.0+(DmInv(2,2,:).*t20)./2.0;
            t2.*t19.*(-1.0./2.0)-(t3.*t18)./2.0;
            t2.*t21.*(-1.0./2.0)-(t3.*t20)./2.0;
            t18.^2./2.0+t20.^2./2.0-1.0./2.0;%constraints
            t19.^2./2.0+t21.^2./2.0-1.0./2.0;
            (t18.*t19)./2.0+(t20.*t21)./2.0];
    
    gradC = pagetranspose([out(1:6,:,:),out(7:12,:,:),out(13:18,:,:)]);
    constraintC = reshape(out(19:end,:,:),3,1,[]);
end

