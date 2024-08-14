function [angle] = dihedralAngle(v1,v2,o1,o2)
%DIHEDRALANGLE compute the dihedralAngle from two triangles separated by an
%edge. v1 & v2 vertices of edge. o1 & o2 vertices opposite of edge 
% 1 cross 1 norm version which should help the symbolic mat lib
    b0 = o1 - v1;
    b1 = v2 - v1;
    b2 = o2 - v2;
    
    b1norm = fastVectorNorm(b1);
    b1normalized = b1./b1norm;

    v = b0 - fastDot(b0,b1normalized)*b1normalized;
    w = b2 - fastDot(b2,b1normalized)*b1normalized;

    x = fastDot(v,w);
    b1nCrossv = fastCross(b1normalized,v);
    y = fastDot(b1nCrossv,w);
    
    angle = atan2(y,x);
end

