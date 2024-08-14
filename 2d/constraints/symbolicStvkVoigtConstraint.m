function [deltax, Minv, gradC, deltaLambda, alphaTilde, constraintC] = symbolicStvkVoigtConstraint(x, DmInv, alphaTilde, lambda, minv)
    minv = 1./mass;
    Minv = diag(minv);

    CgradCSymbols = setupSymbolsGradC(DmInv,x);

    constraintC = CgradCSymbols(19:end);
    gradC = [CgradCSymbols(1:6),CgradCSymbols(7:12),CgradCSymbols(13:18)]';

    topDenominator = gradC(1:2,:)*Minv*gradC(1:2,:)' + alphaTilde(1:2,1:2);
    topNumerator = -constraintC(1:2) - alphaTilde(1:2,1:2)*lambda(1:2)';
    topDeltaLambda = topDenominator\topNumerator;

    denominator = gradC(3,:)*Minv*gradC(3,:)' + alphaTilde(3,3);
    numerator = -constraintC(3) - alphaTilde(3,3)*lambda(3);

    deltaLambda = [topDeltaLambda; numerator/denominator];
    deltax = Minv*gradC'*deltaLambda;
end

function out = setupSymbolsGradC(DmInv_in, x_in)
    x = sym( 'x', [2*3,1], 'real' );
    DmInv = sym( 'DmInv', [2,2], 'real' );
    
    p0 = x(1:2);
    p1 = x(3:4);
    p2 = x(5:6);
    
    Ds = [p0-p2,p1-p2];
    F = Ds*DmInv;
    E = 0.5 *( F' * F - eye(2) );
    constraintC = [E(1,1);E(2,2);E(1,2)];
    
    %need to generate code for all 3 constraints of constraintC
    %here I just show the first one as an example
    gradC = [gradient(constraintC(1),x)';gradient(constraintC(2),x)';gradient(constraintC(3),x)'];

    outsym = [reshape(gradC',[],1);constraintC];
    outs1 = subs(outsym,DmInv, DmInv_in);
    outs2 = subs(outs1,x, x_in);
    out = eval(outs2);
end
