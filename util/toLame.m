function [mu, lambda] = toLame(nu, E)
    % TOLAME Converts Poisson ratio and nu and Young's modulus E to Lamé
    % parameters.
    mu = E / (2 * (1 + nu));
    lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
end

