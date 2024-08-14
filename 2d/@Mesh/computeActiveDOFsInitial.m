function activeDOFs = computeActiveDOFsInitial(mesh)
% active dofs will also have rigid dofs removed (override by adaptive mesh)
    activeDOFs = mesh.unpinnedDOFs;
end