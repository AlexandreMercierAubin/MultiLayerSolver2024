function res = mult3x3TransposeIn3D(A, B)
    %multiply matrix A by B transposed
    %first line
    a = A(1, 1, :) .* B(1, 1, :) + A(1, 2, :) .* B(2, 1, :) + A(1, 3, :) .* B(3, 1, :);
    b = A(1, 1, :) .* B(1, 2, :) + A(1, 2, :) .* B(2, 2, :) + A(1, 3, :) .* B(3, 2, :);
    c = A(1, 1, :) .* B(1, 3, :) + A(1, 2, :) .* B(2, 3, :) + A(1, 3, :) .* B(3, 3, :);
    %second line
    d = A(1, 2, :) .* B(1, 1, :) + A(2, 2, :) .* B(2, 1, :) + A(3, 2, :) .* B(3, 1, :);
    e = A(1, 2, :) .* B(1, 2, :) + A(2, 2, :) .* B(2, 2, :) + A(3, 2, :) .* B(3, 2, :);
    f = A(1, 2, :) .* B(1, 3, :) + A(2, 2, :) .* B(2, 3, :) + A(3, 2, :) .* B(3, 3, :);
    %third line
    g = A(1, 3, :) .* B(1, 1, :) + A(2, 3, :) .* B(2, 1, :) + A(3, 3, :) .* B(3, 1, :);
    h = A(1, 3, :) .* B(1, 2, :) + A(2, 3, :) .* B(2, 2, :) + A(3, 3, :) .* B(3, 2, :);
    i = A(1, 3, :) .* B(1, 3, :) + A(2, 3, :) .* B(2, 3, :) + A(3, 3, :) .* B(3, 3, :);
    
    res(1:3, 1:3, :) = [a, b, c; d, e, f; g, h, i];
end
