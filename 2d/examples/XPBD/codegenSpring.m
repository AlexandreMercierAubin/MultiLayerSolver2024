p1 = sym( 'p1', [2,1], 'real' );
p2 = sym( 'p2', [2,1], 'real' );
k = sym( 'k', 'real' );
l0 = sym( 'l0', 'real' );

norml = sqrt(sum((p2-p1).^2));
C = 0.5*k*(norml-l0)^2;
gradCp1 = gradient(C,p1);
gradCp2 = gradient(C,p2);

%using l to make it easier to read l = norm(p2-p1)
l = sym( 'l', 'real' );
newGradp1 = subs(gradCp1,norml,l)