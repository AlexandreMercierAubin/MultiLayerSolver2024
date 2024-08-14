function isColliding = sphereSphereCollision(spherePos1, spherePos2, radius1, radius2)
%isColliding = sphereSphereCollision(spherePos1, spherePos2, radius1, radius2)    
%assumes the post a column vectors
    dist = spherePos1-spherePos2;
    normDist = fastVecNorm(dist,dist);
    isColliding  = normDist < (radius1+radius2);
end

