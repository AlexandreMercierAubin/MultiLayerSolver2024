function [forces,HessianPsi] = neoGrinspunianBendingEnergy(kd, angleDiff,vertices,normal,normalTilde,area,areaTilde,innerAngle1,innerAngleTilde1,innerAngle2,innerAngleTilde2)
%NEOGRINSPUNIANBENDINGENERGY energy from the new version of bending energy
% it looks awesome compare to the other one............. wooooh!
% Tilde means the second triangle
% https://www-sciencedirect-com.proxy3.library.mcgill.ca/science/article/pii/S1524070313000209
    x1 = vertices(2,:); %middle 1
    x2 = vertices(3,:); %middle 2    
    e0 = (x2-x1);
    norme0 = norm(e0);
    e0normalized = e0./norme0;

    %computing the gradient and hessian
    x0 = vertices(1,:); %opposite 1
    x3 = vertices(4,:);%opposite 2
    e1 = (x0-x2);
    e2 = (x0-x1);
    eTilde1 = (x3-x2);
    eTilde2 = (x3-x1);
    norme2 = norm(e2);
    norme1 = norm(e1);
    
    normeTilde1 = norm(eTilde1);
    normeTilde2 = norm(eTilde2);

    h0 = triangleHeight(area, norme0);
    h1 = triangleHeight(area, norme1);
    h2 = triangleHeight(area, norme2);
    hTilde0 = triangleHeight(areaTilde, norme0);
    hTilde1 = triangleHeight(areaTilde, normeTilde1);
    hTilde2 = triangleHeight(areaTilde, normeTilde2);
    
    omegaTilde = [1/(hTilde0*hTilde0),0,0;
                  1/(hTilde1*hTilde0),1/(hTilde1*hTilde1),0;
                  1/(hTilde2*hTilde0),1/(hTilde2*hTilde1),1/(hTilde2*hTilde2)];
    omegaTilde = omegaTilde + triu(omegaTilde.',1);%symetric copy of the top

    omega = [1/(h0*h0),0,0;
             1/(h1*h0),1/(h1*h1),0;
             1/(h2*h0),1/(h2*h1),1/(h2*h2)];
    omega = omega + triu(omega.',1);
  
    e1normalized = e1./norme1;
    e2normalized = e2./norme2;
    eTilde1normalized = eTilde1./normeTilde1;
    eTilde2normalized = eTilde2./normeTilde2;

    %oof, I had to visualize them all with the right hand rule to make sure
    m0 = cross(e0normalized,normal);
    m1 = cross(e1normalized,normal);
    m2 = cross(normal,e2normalized);
    
    mTilde0 = cross(normalTilde,e0normalized);
    mTilde1 = cross(normalTilde,eTilde1normalized);
    mTilde2 = cross(eTilde2normalized,normalTilde);

    psiPrime = -2*kd*angleDiff;
    psiPrimePrime = -2*kd; %Flipping the hessian makes this work for some reason

    gradThetax0 = -(1/h0)*normal';
    gradThetax1 = (cos(innerAngle2)/h1)*normal' + (cos(innerAngleTilde2)/hTilde1)*normalTilde';
    gradThetax2 = (cos(innerAngle1)/h2)*normal' + (cos(innerAngleTilde1)/hTilde2)*normalTilde';
    gradThetax3 = -(1/hTilde0)*normalTilde';

    gradTheta = [gradThetax0;gradThetax1;gradThetax2;gradThetax3]';

    forces = psiPrime * gradTheta';

    HessOuterProduct = gradTheta'*gradTheta;

    HessianAngleComponent = zeros(12,12);

    id0 = 1:3;
    id1 = 4:6;
    id2 = 7:9;
    id3 = 10:12;

    %H00
    M0 = (normal'*m0);
    M1 = (normal'*m1);
    M2 = (normal'*m2);
    Q0 = omega(1,1)*M0; %carefull, omega indexing starts at 1 because of matlab, but paper says zero
    HessianAngleComponent(id0,id0) = -S(Q0);

    %H33
    MTilde0 = (normalTilde'*mTilde0);
    MTilde1 = (normalTilde'*mTilde1);
    MTilde2 = (normalTilde'*mTilde2);
    QTilde0 = omegaTilde(1,1)*MTilde0; 
    HessianAngleComponent(id3,id3) = -S(QTilde0);

    %H11
    squared0norm = e0*e0'; %squared norm
    N0 = M0./squared0norm;
    P11 = omega(2,2)*cos(innerAngle1)*M1';
    NTilde0 = MTilde0./squared0norm;
    PTilde11 = omegaTilde(2,2)*cos(innerAngleTilde1)*MTilde1';
    N0N0Tilde = N0 + NTilde0;
    HessianAngleComponent(id1,id1) = S(P11) + S(PTilde11) - N0N0Tilde;

    
    %H22
    P22 = omega(3,3)*cos(innerAngle2)*M2';
    PTilde22 = omegaTilde(3,3)*cos(innerAngleTilde2)*MTilde2';
    HessianAngleComponent(id2,id2) = S(P22) + S(PTilde22) - N0N0Tilde;

    %H10
    P10 = omega(2,1)*cos(innerAngle1)*M0';
    Q1 = omega(1,2)*M1;
    HessianAngleComponent(id1,id0) = P10 - Q1;

    %H20
    P20 = omega(3,1)*cos(innerAngle2)*M0';
    Q2 = omega(1,3)*M2;
    HessianAngleComponent(id2,id0) = P20 - Q2;

    %H13'
    PTilde10 = omegaTilde(2,1)*cos(innerAngleTilde1)*MTilde0';
    QTilde1 = omegaTilde(1,2)*MTilde1;
    H13 = PTilde10 - QTilde1;
    HessianAngleComponent(id3,id1) = H13';

    %H23'
    PTilde20 = omegaTilde(3,1)*cos(innerAngleTilde2)*MTilde0';
    QTilde2 = omegaTilde(1,3)*MTilde2;
    H23 = PTilde20 - QTilde2;
    HessianAngleComponent(id3,id2) = H23';

    %H12'
    P12 = omega(2,3)*cos(innerAngle1)*M2';
    P21 = omega(3,2)*cos(innerAngle2)*M1';
    PTilde12 = omegaTilde(2,3)*cos(innerAngleTilde1)*MTilde2';
    PTilde21 = omegaTilde(3,2)*cos(innerAngleTilde2)*MTilde1';
    H12 = P12 + P21' + PTilde12 + PTilde21' + N0N0Tilde;
    HessianAngleComponent(id2,id1) = H12';

    HessianAngleComponentSymetric = triu(HessianAngleComponent.',1) + tril(HessianAngleComponent);
    HessianPsi = (psiPrime*HessianAngleComponentSymetric + psiPrimePrime*HessOuterProduct);

    function [angleTri] = angleTriangleA(edgea,edgeb,edgec)
        numerator = edgeb*edgeb + edgec*edgec - edgea*edgea;
        denominator = 2*edgeb*edgec;
        angleTri  = acos(numerator/denominator);
    end    
    
    function Smatrix = S(A)
        Smatrix = A + A';
    end

    function area = triangleArea(edgeLength1, edgeLength2, angleInBetween)
        area = 0.5 * edgeLength1 * edgeLength2 * sin(angleInBetween);
    end    
    
    function curvature = locallyIntegratedMeanCurvature(angle)
        curvature = 2*tan(angle/2);
    end

    function height = triangleHeight(area, edgeLength)
        height = 2 * area / edgeLength;
    end

end

