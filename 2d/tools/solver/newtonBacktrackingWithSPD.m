function [x, it] = newtonBacktrackingWithSPD(f, x, maxIterations, gamma, c, eps)
    % notation from https://sites.math.washington.edu/~burke/crs/516/notes/backtracking.pdf

    if nargin < 3
        maxIterations = 10;
    end
    
    if nargin < 4
        gamma = 0.5;
    end
    
    if nargin < 5
        c = 1e-4;
    end
    
    if nargin < 6
        eps = 1e-9;
    end
    
    fx = f(x);
    normFx = norm(fx);
    
    for it = 1:maxIterations
        J = zeros(numel(fx), numel(x));
        for i = 1:numel(x)
            delta = zeros(numel(x), 1);
            delta(i, 1) = eps;
            J(:, i) = (f(x + delta) - fx) / eps;
        end
        
        J = PositiveSemiDefiniteProjection(J);
        
        d = -(J \ fx);
      
        t = 1;
        
        while normFx - norm(f(x + t * d)) < c * t * normFx
            t = gamma * t;
        end
        
        x = x + t * d;
        
        fx = f(x);
        normFx = norm(fx);
        if normFx < eps
            return;
        end
    end
end

