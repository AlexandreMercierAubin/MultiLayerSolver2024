function n_vec = computeFnormalComponent(meshes,F,settings)
    numt = size(meshes.t,1);
    if isa(meshes,'AdaptiveMesh3D') 
        numt = numel(meshes.ElasticTetInds);
    end    
    n_vec = zeros(9*numt,1);
    if size(meshes.t,2) == 3 && settings.addShellNormalDeformation
        % we need to add the normal change to the deformation gradient of shells
        position = meshes.B'*F;
        formattedPositions = reshape(position',3,[])';

        n = normals(formattedPositions, meshes.t);
        triNormals = normr(n);
        normalStack = reshape(triNormals',[],1);
        n_vec = meshes.referenceSpaceNormalAssembled*normalStack;
    end
end