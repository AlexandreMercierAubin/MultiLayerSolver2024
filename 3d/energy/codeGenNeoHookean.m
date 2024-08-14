clear;
dimension = 3;
F = sym( 'F', [dimension,dimension], 'real' );
mu = sym( 'mu', 'real'  );
lambda = sym( 'lambda', 'real'  );
volume = sym( 'volume', 'real' );

J = det(F);
C = F'*F;
% Cbar = (1/J) * C;
Ic = trace(C);
alpha = 1 + mu/lambda - mu/(4*lambda);
Ja = J - alpha;

%https://graphics.pixar.com/library/StableElasticity/paper.pdf eq 14
psi = 0.5 * mu * (Ic-dimension) + 0.5 * lambda *Ja*Ja - 0.5 * mu * log(Ic + 1); 

% psi = 0.5 * mu * (trace(Cbar) - dimension) + lambda * (J-1)^2;
assume(psi,'real');

psi = psi*volume;
matlabFunction( psi, 'File', 'NeoHookeanPsi.c' );

dpsidF = -gradient(psi,F(:));
d2psidF2 = -hessian(psi,F(:));

d2psidF2 = d2psidF2(:);

% asking matlab to make a file for all parts simultaneously allows for a
% better use (elimination) of common sub expressions.  Here, the force is 
% put following Hessian so it is easy to pull off the end of the code gen.
out = simplify([ reshape(d2psidF2,[],1); dpsidF ; psi],'Criterion', 'preferReal'); 
ccode( out, 'File', 'NeoHookeanPsiGradHess.c' );