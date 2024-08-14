dimension = 3;
F = sym( 'F', [dimension,dimension], 'real' );
materialMu = sym( 'materialMu', 'real'  );
materialLambda = sym( 'materialLambda', 'real'  );
lambda = sym( 'lambda', 'real'  );
volume = sym( 'volume', 'real' );
alphaTilde = sym( 'alphaTilde', 'real' );
gamma = sym( 'gamma', 'real' );

E = 0.5 *( F' * F - eye(dimension) );
froE2 = sum( E.*E, 'all');
trE = trace(E);
psi = materialMu * froE2 + materialLambda * 0.5 * trE*trE;
constraintC = mesh3d.elV(tet) * psi;

gradC = gradient(constraintC);
deltaLambda = (-constraintC - alphaTilde*lambda - gamma*constraintC)/()