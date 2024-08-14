function computeRigidForce(mesh, h)
    for i = 1:1:numel(mesh.RigidBodies)
        mesh.RigidBodies(i).computeRigidForce(h);
    end
end
