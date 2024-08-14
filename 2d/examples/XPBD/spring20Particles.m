m = 1;
x = zeros(20,2);
v = zeros(20,2);
springs = [1:19;2:20]';

g= -9.8;
f = [0,g/m];
alpha = 10^(-8);

x(:,1) = 0:1:19;
h = 1/60 ;

for i = 1:200
 [x,v] = step(h,x,v,springs,f,alpha,m);
 plot(x(:,1),x(:,2),'b-o');
 axis([-20,20,-20,2])
 pause(1/60);
end
tmpEnd = 1;

function [x,v] = step(h,xold, vold, springs, f , alpha,m)
lambda = zeros(19,1);
v = vold + h*f;
x = xold + h*v;
x(1,:) = [0,0];

for it = 1 : 100
    for s = 1:size(springs,1)
        p1 = springs(s,1);
        p2 = springs(s,2);
        [deltaLambda, dx1, dx2,C] = distanceConstraint.projectNoBeta(h, x(p1,:), x(p2,:),m,m,lambda(s),alpha,1);
        x(p1,:) = x(p1,:) + dx1;
        x(p2,:) = x(p2,:) + dx2;
        lambda(s) = lambda(s) + deltaLambda;
    end

    x(1,:) = [0,0];
end

v = (x - xold)./h;

end