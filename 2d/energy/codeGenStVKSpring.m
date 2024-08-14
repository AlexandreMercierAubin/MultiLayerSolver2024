clear;
x = sym( 'x', [4,1], 'real' );
mass = sym( 'mass', [4,1], 'real' );
k = sym( 'k' , 'real' );
l0 = sym( 'l0' , 'real' );
alpha = sym( 'alpha' , 'real' );
h = sym( 'h' , 'real' );
lambda = sym( 'lambda' , 'real' );

p0 = x(1:2);
p1 = x(3:4);

direction = p1-p0;
l = sqrt(sum((direction).^2));
constraintC = 0.5*k*(l-l0); %xpbdParallel eq 17


gradC = gradient(constraintC,x); %-------
gradp0 = -0.5*k*direction/(l);