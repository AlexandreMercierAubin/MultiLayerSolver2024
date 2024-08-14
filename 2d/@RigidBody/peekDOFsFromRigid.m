function [p,v] = peekDOFsFromRigid(body, Angle, Position, h)
    %SETDOFSFROMRIGID peeks at the mesh particles that are part of this body
    %using this body's position and orientation without updating the rigid
    %body
    rot = [cos(Angle), -sin(Angle); sin(Angle), cos(Angle)];

    inds = body.Indices;
    inds2 = inds*2;
    inds1 = inds2 - 1;
    
    pos = Position;
    vertexDisp = body.Mesh.VertexDisp(inds, :);
    Rpx = vertexDisp(:, :) * rot(1, :)';
    Rpy = vertexDisp(:, :) * rot(2, :)';

    p = body.Mesh.p;
    prevp = [body.Mesh.p(inds1),body.Mesh.p(inds2)]';
    p(inds1) = pos(1) + Rpx;
    p(inds2) = pos(2) + Rpy;
    
    v = body.Mesh.v;
    vel = ([body.Mesh.p(inds1),body.Mesh.p(inds2)]' - prevp)/h;
    v(inds1) = vel(1,:);
    v(inds2) = vel(2,:);
end