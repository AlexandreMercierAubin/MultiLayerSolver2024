function x = solveBGS(iterations, block_size, A, b, initial_x)
    n = size(b, 1);
    
    if nargin > 4
        x = initial_x;
    else
        x = zeros(n, 1);
    end
    
    for it = 1:iterations
        for i = 1:block_size:n
            ix = i:i + block_size - 1;
            x(ix) = inv(A(ix, ix)) * (b(ix) - A(ix, :) * x + A(ix, ix) * x(ix));
        end
    end
end

