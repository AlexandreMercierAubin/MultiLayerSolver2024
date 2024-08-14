clear;
dimension = 3;
F = sym( 'F', [dimension,dimension], 'real' );
mu = sym( 'mu' , 'real' );
lambda = sym( 'lambda' , 'real' );
volume = sym( 'volume' , 'real' );

E = 0.5 *( F' * F - eye(dimension) );
froE2 = trace(E*E');
trE = trace(E);
psi = mu * froE2 + lambda * 0.5 * trE*trE;
psi = psi*volume;
matlabFunction( psi, 'File', 'STVK.c' );

dpsidF = -gradient(psi,F(:));

d2psidF2 = -hessian(psi,F(:));
d2psidF2 = d2psidF2(:);

% asking matlab to make a file for all parts simultaneously allows for a
% better use (elimination) of common sub expressions.  Here, the force is 
% put following Hessian so it is easy to pull off the end of the code gen.
out = simplify([ reshape(d2psidF2,[],1); dpsidF ; psi],'Criterion', 'preferReal'); 
ccode( out, 'File', 'STVKGradHess.c' );