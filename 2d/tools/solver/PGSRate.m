function maxEig = PGSRate(M, inds)
    n = size(M, 1);
    PGS = zeros(n, n);
    for i = 1:n
        x = zeros(n, 1);
        x(i) = 1;
        PGS(:, i) = oneStepPGS(M, x, inds);
    end
    
    maxEig = max(abs(eig(PGS)));
end


function xkp1 = oneStepPGS(A, xk, inds)
    n = size(A, 1);
    xkp1 = xk;
    for i = 1:n
        xkp1(i) = (-A(i, :) * xkp1 + A(i, i) * xkp1(i)) / A(i, i);
    end
    
    z = zeros(n, 1);
    xkp1(inds) = max(z(inds), xkp1(inds));
end