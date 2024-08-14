function [deltaLambda, deltaxElastic, rigidCOM]=rigidBoundary1D(h, pElastic, pRigid, massElastic, lambda, alpha, rigidBodyMass, iterateRigidCOM)
    assert(size(pRigid,1)==1)
    assert(size(pElastic,1)==1)
    %for the mass splitting term nb, we divide the mass of the
    %constraint by the number of constraints on the rigid body, but
    %let the impulse use the full mass.
    minv = 1/massElastic;
    %this is 1D so inertia and nb are always 1
    iterateInertia = 1;
    nb = 1;

    direction = pElastic-pRigid;
    constraintC = norm(direction);
    if constraintC == 0
        %no need to solve, the constraint is exactly satisfied
        deltaLambda = 0;
        deltaxElastic = zeros(size(pElastic));
        rigidCOM = iterateRigidCOM;
        return;
    end
    n = direction./constraintC;

    Iinv = 1/iterateInertia;
    invMassRigid = 1/rigidBodyMass;
    w = invMassRigid/nb; %it says split mass and inertia in the mass splitting paper

    minvSum = minv+w;
    
    h2 = h*h;
    alphaTilde = alpha/h2;

    numerator = -constraintC - alphaTilde*lambda;
    denumerator = minvSum + alphaTilde;
    deltaLambda = numerator/denumerator;

    p = deltaLambda*n;
    deltaxElastic = minv*p;
    deltaCOM = (p.*invMassRigid);
    rigidCOM = iterateRigidCOM - deltaCOM;
end