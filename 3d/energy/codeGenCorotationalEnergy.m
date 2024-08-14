clear;
dimension = 3;
F = sym( 'F', [dimension,dimension], 'real' );
mu = sym( 'mu', 'real' );
lambda = sym( 'lambda', 'real' );
volume = sym( 'volume', 'real' );

%cheap S of svd decomposition
S = eig(F'*F);
sumS = sum((S-ones(dimension,1)).^2);
prodS = sum(S)-dimension;
psi = mu *sumS + 0.5*lambda*prodS*prodS;

psi = psi*volume;
assume(psi, 'real');

dpsidF = -gradient(psi,F(:));

d2psidF2 = -hessian(psi,F(:));
d2psidF2 = d2psidF2(:);

% asking matlab to make a file for all parts simultaneously allows for a
% better use (elimination) of common sub expressions.  Here, the force is 
% put following Hessian so it is easy to pull off the end of the code gen.
% out = simplify([ reshape(d2psidF2,[],1); dpsidF ; psi],'Criterion', 'preferReal'); 
out = [ reshape(d2psidF2,[],1); dpsidF ; psi]; 


ccode( out, 'File', 'corotational.c' ); 
% save('corot.c', 'out');
