function [E,nu] = toEmu(lambda,mu)
    %TOEMU converts lamé parameters to young modulus E and poison ratio nu
    E =  mu*(3*lambda + 2*mu)/(lambda + mu);
    nu = lambda/(2*(lambda+mu));
end

