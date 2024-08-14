function x = solveGS(iterations, A, b, initial_x)
    n = size(b, 1);
    
    if nargin > 3
        x = initial_x;
    else
        x = zeros(n, 1);
    end
    
    for it = 1:iterations
        for i = 1:n
            sigma = dot(A(i, 1:i - 1)', x(1:i - 1)) + dot(A(i, i + 1:n)', x(i + 1:n));
            x(i) = (b(i) - sigma) / A(i, i);
        end
    end
end

