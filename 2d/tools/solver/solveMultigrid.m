function [x, resInterpolated] = solveMultigrid(A, x0, b, presmoothing, postsmoothing)
%[res,x] = solveMultigrid(A,x0,b)
% res: residual error
% x: solution
% https://tiantianliu.cn/papers/xian2019multigrid/xian2019multigrid.html

%     smoother = @jacobi;
%     for i = 1:presmoothing
%         u = smoother(A, x, b, res);
%     end
    
    %coarse grid correction
    x = x0;
    n = numel(x0);

    if n > 2 && mod(n,2)==0
        %compute residual
        res = b-A*x;
        
        %restrict
        resRestricted = res(1:2:end);
        ARestricted = A(1:2:end,1:2:end);
        bRestricted = b(1:2:end);
        
        
        [xTmp, resCoarse] = solveMultigrid(ARestricted, x(1:2:end), bRestricted, presmoothing, postsmoothing);

        %interpolate
        resInterpolated = zeros(numel(resCoarse)*2,1);
        resInterpolated(1:2:end) = resCoarse;
        resInterpolated(2:2:end-1) = 0.5*(resInterpolated(1:2:end-2)+resInterpolated(3:2:end));
        
        %improve solution
        x = x + resInterpolated;
    else
        x = A\b;
        resInterpolated = b-A*x;
    end
    
%     for i = 1:postsmoothing
%         u = smoother(A, x, b, res);
%     end

end

