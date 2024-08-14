function setDOFsFromRigid( body, h )
    %SETDOFSFROMRIGID updates the mesh particles that are part of this body
    %using this body's position and orientation
    [pList,vList,inds] = peekDOFsFromRigid( body, h, body.Rotation, body.Position, body.Mesh.p);
    body.Mesh.p(inds) = pList;
    body.Mesh.v(inds) = vList;
end