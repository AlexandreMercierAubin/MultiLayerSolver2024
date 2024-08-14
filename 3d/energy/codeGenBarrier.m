clear;
s = 3;
x = sym( 'x', [s,1], 'real' );
c = sym( 'c', [s,1], 'real' );
dhat = sym('dhat', 'real');

d=sqrt(sum((x-c).^2));
B = -(d-dhat)^2*log(d/dhat);

gradientB = simplify(gradient(B,x));
hessianB = jacobian(gradientB,x);

out = [ reshape(hessianBt,[],1); gradientBt]; 
ccode( out, 'File', 'IPC.c' );

