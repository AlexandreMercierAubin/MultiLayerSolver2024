function [angle] = symDihedralAngleAtan(v1,v2,o1,o2)
%DIHEDRALANGLE compute the dihedralAngle from two triangles separated by an
%edge. v1 & v2 vertices of edge. o1 & o2 vertices opposite of edge 
    b1 = o1 - v1;
    b2 = v2 - v1;
    b3 = o2 - v1;
    
    n1 = fastCross(b1,b2);
    n1norm = fastVectorNorm(n1);
    n1 = n1./n1norm;
    
    n2 = fastCross(b2,b3);
    n2norm = fastVectorNorm(n2);
    n2 = n2./n2norm;
    
    b2norm = fastVectorNorm(b2);
    e = b2./b2norm;
    
    angle = dihedralAngleFromNormals(n1',n2',e');
    %for some reason atan2 does not work well todo: code atan here and try
    %to  bake atan2 in the mex function?
    % https://www.mathworks.com/help/symbolic/atan2.html#:~:text=atan2(%20Y%20%2C%20X%20)%20computes,computes%20arctangents%20element%20by%20element.
end

