function x = truncSVDSolve(A, b, threshold)
    [U, S, V] = svd(A);

    sd = diag(S);

    sd(sd < threshold) = Inf;

    invS = diag(1 ./ sd);

    Ainv = V * invS * U';

    x = Ainv * b;
end