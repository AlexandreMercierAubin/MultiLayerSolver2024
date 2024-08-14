clear;
s = 4;
mu = sym( 'mu' , 'real');
lambda = sym( 'lambda' , 'real');
gforce = sym( 'gforce', 'real');
d_hat = sym( 'd_hat' , 'real');
C_hat = sym( 'C_hat' , 'real');
e = 1;
x = sym( 'x', [s,1], 'real' );
v = sym( 'v', [s,1], 'real' );
b = sym('b', 'real');
h = sym('h', 'real');
K = sym('K', 'real');
M = sym('M', [s,s], 'real');
Minv = sym('Minv', [s,s], 'real');

fd = sym( 'fd', [s,1], 'real' );
fe = sym( 'fe', [s,1], 'real' );

B= sym( 'B', [9*e,s], 'real' );
Bx = B*x;
for i = 1:e
    F = reshape(Bx(9*i-8:9*i),3,3);
    

    J = det(F);
    C = F'*F;
    Ic = trace(C);

    %https://graphics.pixar.com/library/StableElasticity/paper.pdf eq 14
    Psi(i) = 0.5 * mu * (Ic-3) + 0.5 * lambda *(J-1)*(J-1);
end
fe(3) = gforce;

x_hat = x + h*v + h*h*Minv*fe;
xdiff = (x-x_hat);

%h*h*psi, but we are using elastic forces here instead so h*psi
E = 0.5*xdiff'*M*xdiff - h * h * x' * fd + h *h * sum(Psi);

B_t = E + K * sum(b);

gradientBt = gradient(B_t,x);
hessianBt = jacobian(gradientBt,x);

out = [ reshape(hessianBt,[],1); gradientBt]; 
ccode( out, 'File', 'IPC.c' );

