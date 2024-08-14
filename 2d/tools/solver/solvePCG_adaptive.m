function [x,resvec,xarray] = solvePCG( iterations, A, b, initial_x, adaptive_preconditioner )
% SOLVEPCG Conjugate Gradient on Ax=b, using a preconditioner and starting 
% with the initial point (should be zero! Use shifted Krylov to warm start
% instead).  
%   iterations  maximum number of iteraitons
%   A           matlab function that multiplies by A
%   b           the right hand side
%   initial_x   should be zero! 
%   cache       has the preconditioner, a matlab function that applies the 
%               preconditioner. 
    x = initial_x;
    r = b - A(x);
    z = adaptive_preconditioner.precondition(r);
    d = z;
    rTz = r' * z;

    if nargout > 1
        resvec = zeros(iterations,1);
        if nargout > 2
            xarray = zeros(numel(x),iterations);
        end
    end

    for i=1:iterations
        adaptive_preconditioner.checkUpdate(i);
        Ap = A(d);
        qAp = (d'*Ap);
        if qAp == 0
            alpha = 0;
        else
            alpha = rTz / (d'*Ap);
        end
        x = x + alpha * d;
        r = r - alpha * Ap;
%         if norm(r) < 1e-30
%               break;
%         end
        z = preconditioner(r);
        rTznew = r'*z;
        beta = rTznew / rTz;
        rTz = rTznew;
        d = z + beta * d;
        if nargout > 1
            resvec(i) = norm(A(x)-b);
            if nargout > 2
                xarray(:,i)= x;
            end
        end
    end
end