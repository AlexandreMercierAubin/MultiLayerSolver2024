function vec = diagonalDominanceVector(A)
    vec = zeros(size(A, 1), 1);
    for i = 1:size(A, 1)
        vec(i) = abs(A(i, i)) - sum(abs(A(i, [1:i - 1, i + 1:end])));
    end
end