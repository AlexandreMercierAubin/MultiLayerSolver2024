function Ahat = nearestSPDS(A)
% finds an spd projection of a sparse matrix
% arguments: (input)
%  A - square matrix, which will be converted to a Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.

if nargin ~= 1
  error('Exactly one argument must be provided.')
end

sigma = 0;
Ahat = A;
Ahat = (Ahat + Ahat')/2;

% issymmetric(Ahat)
[u,v] = eigs(Ahat,1,"smallestreal");

% while(v < sigma)
%     Ahat = Ahat - u*u'*v;
%     Ahat = (Ahat + Ahat')/2;
%     [u,v] = eigs(Ahat,1,"smallestreal");
% end

% clamp the eigs that haven't converged yet.
[V,D] = eig(full(Ahat));
u = diag(D);
u(u < sigma) = sigma;
Ahat = V*diag(u)*V';

Ahat = sparse(Ahat);






