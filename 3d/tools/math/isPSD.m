function [isMatrixPSD, isMatrixSemiDef] = isPSD(A, tol)
    if nargin < 2
        tol = 1e-8;
    end
    tf = issymmetric(A);
    d = eig(A);
    isposdef = all(d > tol);
    isMatrixPSD = tf && isposdef;
    issemidef = all(d > -tol);
    isMatrixSemiDef = tf && issemidef;
end

