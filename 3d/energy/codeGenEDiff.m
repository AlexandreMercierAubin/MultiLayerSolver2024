clear;
dimension = 3;
F = sym( 'Fa', [dimension,dimension], 'real' );
oldF = sym( 'Fb', [dimension,dimension], 'real' );
h = sym( 'h' , 'real' );

E = 0.5 *( F' * F - eye(dimension) );
oldE = 0.5 *( oldF' * oldF - eye(dimension) );
EDiff = (E - oldE)/h;

EdiffSquared = EDiff.*EDiff;

ccode( sum(EdiffSquared(:)), 'File', 'EDiff.c' );