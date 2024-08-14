function f = peekCurrentForce(mesh, h, AngularVelocity, RigidBodyPosition, Inertia, meshf, meshp)
    rigid = zeros(0,0);
    Force = zeros(3,numel(mesh.RigidBodies));
    Torque = zeros(3,numel(mesh.RigidBodies));
    for i = 1:1:numel(mesh.RigidBodies)
        [force,torque] = mesh.RigidBodies(i).peekRigidForce(h, AngularVelocity(:,i), RigidBodyPosition(:,i), Inertia(:,:,i), meshf, meshp);
        Force(:,i) = force;
        Torque(:,i) = torque;
    end
    if numel(mesh.RigidBodies) > 0
        rigid = zeros(6*numel(mesh.RigidBodies),1);
        rigid(4:6:end,1) = Torque(1,:);
        rigid(5:6:end,1) = Torque(2,:);
        rigid(6:6:end,1) = Torque(3,:);
        rigid(1:6:end,1) = Force(1,:);
        rigid(2:6:end,1) = Force(2,:);
        rigid(3:6:end,1) = Force(3,:);
    end
    f = [meshf(mesh.ElasticDOFs); rigid];
    
end
