function newxi = rotateAllVerts2D(rigidIDbyVert, ri, R, COM, xi)
    rigidVerts = find(rigidIDbyVert ~= 0);
    riCropped = permute(ri(rigidVerts,:),[3,2,1]);
    queriedxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,rigidIDbyVert(rigidVerts)), COM(:,:,rigidIDbyVert(rigidVerts)));
    newxi = xi;
    newxi([rigidVerts*2-1,rigidVerts*2]') = queriedxi;
end