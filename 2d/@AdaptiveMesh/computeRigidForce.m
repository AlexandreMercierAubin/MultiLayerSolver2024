function computeRigidForce(mesh,h)
    %for loop that call the internal computeRigidForce of each rigidbodies
    for i = 1:1:size(mesh.RigidBodies,2)
        mesh.RigidBodies(i).computeRigidForce(h);
    end
end
