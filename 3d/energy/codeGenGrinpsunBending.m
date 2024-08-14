clear;
vertices = sym( 'vertices', [12,1], 'real' );
v1 = vertices(1:3);
v2 = vertices(4:6);
o1 = vertices(7:9);
o2 = vertices(10:12);
angleRest = sym( 'angleRest', 'real' );
kd = sym( 'kd', {'positive','real'} );

%We have to manually add the check for the normal and rest angle inversion when the input to acos is <-1 || >1
angle = symDihedralAngleAtan(v1,v2,o1,o2); 

assume(angle , 'real');
angleDiff = (angle-angleRest);
psi = kd*angleDiff*angleDiff;
assume(psi , 'real');
% psi = simplify(expand(psi));

gradientEdge = gradient(psi, vertices);

hessianEdge = jacobian(gradientEdge, vertices);

out1 = [reshape(hessianEdge,[],1); gradientEdge; angle; psi ];
assume(out1 , 'real');
% simplify(out1,1);
matlabFunction(out1, 'File', 'bendingEnergy.m','Optimize',true)
% ccode( out1, 'File', 'bendingAtanv3.c');