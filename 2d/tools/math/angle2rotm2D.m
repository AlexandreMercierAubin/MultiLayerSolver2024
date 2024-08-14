function R = angle2rotm2D(angle)
%makes a rotation matrix from an angle in radians
    angle=reshape(angle,1,1,[]);
    R = [cos(angle),-sin(angle);sin(angle),cos(angle)];
end
