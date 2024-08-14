function [n1] = triNormal(v1,v2,v3)
%finds the normal of a triangle according to its 3 vertices
    b1 = v3 - v1;
    b2 = v2 - v1;
    
    n1 = fastCross(b1,b2);
    n1norm = fastVectorNorm(n1);
    n1 = n1./n1norm;
end

