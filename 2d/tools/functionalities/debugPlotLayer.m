function debugPlotLayer(mesh2d,layer)
    pr = reshape(mesh2d.p, 2, mesh2d.N)' + mesh2d.plotOffset;
    tintCol = [0.8 0.5 0.2];
    defaultBlue = [0.3,0.46,0.8];
    mesh2d.RenderPatch.Vertices = pr;
    isRigid = layer~=0;

    isRigidVert = false(mesh2d.N,1);
    isRigidVert(mesh2d.t(isRigid,:)) =  true;
    implicitelyRigidElements = all(isRigidVert(mesh2d.t));
    isRigidTotal = implicitelyRigidElements | isRigid;

    mesh2d.RenderPatch.FaceVertexCData(isRigidTotal,:)=repmat(tintCol,sum(isRigidTotal),1);
    mesh2d.RenderPatch.FaceVertexCData(layer==0,:)=repmat(defaultBlue,sum(~isRigidTotal),1);
end

