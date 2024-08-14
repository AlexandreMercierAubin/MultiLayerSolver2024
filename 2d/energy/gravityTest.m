xPrev = sym( 'xPrev', [2*3,1], 'real' );
x = sym( 'x', [2*3,1], 'real' );
mass = sym( 'mass', [2*3,1], 'real' );
g = sym( 'g', 'real' );
h = sym( 'h', 'real' );

deltax = x-xPrev;

constraintC = sum((deltax(2:2:end) - mass(2:2:end)*g*h^2).^2);
gradientC = gradient(constraintC,x);

constraintC2 = mass(2:2:end)'*g*h*deltax(2:2:end);
gradientC2 = gradient(constraintC2,x);