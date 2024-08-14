function [isColliding, distCenter, centers, radius] = AABBcubeCollision(pos1min, pos1max,pos2min,pos2max)
  % this detects collision of two axis aligned boxes defined by their mins
  % and maxs, finds the center and adjust the AABB so it forms a cube
  % instead of an arbitrary box
%     vecRes = (pos1min <= pos2max & pos1max >= pos2min);
    c1 = (pos1min+pos1max)./2;
    c2 = (pos2min+pos2max)./2;
    r1 = c1-pos1min;
    r2 = c2-pos2min;
    maxr1 = max(r1);
    maxr2 = max(r2);
    r1(:) = maxr1;
    r2(:) = maxr2;

    collisionChecks = abs(c1-c2) > (r1+r2);
    if nargout > 1
        distCenter = norm(c1-c2);
        centers = [c1;c2];
        radius = [r1;r2];
    end
    
    isColliding  = all(~collisionChecks);
end

