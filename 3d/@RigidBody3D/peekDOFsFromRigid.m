function [pList,vList,inds] = peekDOFsFromRigid( body, h, R, bodyp, mesh3Dp)
    %SETDOFSFROMRIGID updates the mesh particles that are part of this body
    %using this body's position and orientation
    
    inds3 = body.Indices * 3;
    inds2 = inds3-1;
    inds1 = inds3-2;
    inds = reshape([inds1;inds2;inds3],[],1);

    prevP = [mesh3Dp(inds1),mesh3Dp(inds2),mesh3Dp(inds3)]';
    newP = R * body.Mesh.VertexDisp(body.Indices, :)' + bodyp;
    pList = reshape(newP,[],1);
    
    newP = body.Rotation * body.Mesh.VertexDisp(body.Indices, :)' + body.Position;

    vel = (newP-prevP)/h;
    vList = reshape(vel,[],1);
end
