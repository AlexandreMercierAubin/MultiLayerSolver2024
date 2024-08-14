function [isColliding,distCenter, centers, radius] = AABBboxesCollision(pos1min, pos1max,pos2min,pos2max)
  % this detects collision of two axis aligned boxes defined by their mins
  % and maxs
%     vecRes = (pos1min <= pos2max & pos1max >= pos2min);
    c1 = (pos1min+pos1max)./2;
    c2 = (pos2min+pos2max)./2;
    r1 = c1-pos1min;
    r2 = c2-pos2min;
    collisionChecks = abs(c1-c2) > (r1+r2);
    if nargout > 1
        distCenter = norm(c1-c2);
        centers = [c1;c2];
        radius = [r1;r2];
    end

    isColliding  = all(~collisionChecks);
end

