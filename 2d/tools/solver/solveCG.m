function [x,resvec,xarray, d] = solveCG(iterations, A, b, initial_x, d_in)
    x = initial_x;
    r = b - A * x;
    if nargin <5 || isempty(d_in)
        d = r;
    else
        d = d_in;
    end
    
    if nargout > 1
        resvec = zeros(iterations,1);
        if nargout > 2
            xarray = zeros(numel(x),iterations);
        end
    end

    for i = 1:iterations
        alpha = (r' * r) / (d' * A * d);
        x = x + alpha * d;
        
        oldr = r;
        r = r - alpha * A * d;
        beta = (r' * r) / (oldr' * oldr);
        d = r + beta * d;
        if nargout > 1
            resvec(i) = norm(A*x-b);
            if nargout > 2
                xarray(:,i)= x;
            end
        end
    end
end

