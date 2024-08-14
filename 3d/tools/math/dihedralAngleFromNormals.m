function [angle] = dihedralAngleFromNormals(N1,N2,EV)
%DIHEDRALANGLE compute the dihedralAngle from two normals separated by an
%edge. 
    NCross = fastCross(EV,N1);
    NCrossNorm = NCross./fastVectorNorm(NCross);
    E = [EV', N1', NCrossNorm];
    c = E'*N2';
    angle = atan2(c(3),c(2));
%     angle = acos(dot(N1,N2));
end

