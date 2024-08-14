function setRigidMotion( mesh, angularVelocity, velocity, h )
%SETRIGIDMOTION Sets mesh point velocities to have given velocity
%   mesh.rigidTransform( omega, velocity )
%   Angular motion is about the center of mass
%   Both omega and v and 3 component vectors in the global frame

    [R, COMinitial, inertia, ~, ri, ~, ~, ~]=AdaptiveMesh3D.makeRigidsFromConnectivity(1, mesh.p, zeros(size(mesh.p)), ones(mesh.N,1), false, ones(mesh.N,1), ones(size(mesh.t,1),1),[], mesh.mass, true(size(mesh.p)));
    [R, COM, ~] = RigidBody3D.rigidSymplecticImperativeStep(COMinitial, R, h, velocity, angularVelocity, inertia);
    % xi = xpbdLayer3D.rotateAllVerts(ones(mesh.N,1), ri, R, COM, mesh.p);

    linearDx = h*velocity;
    mesh.prevp = xpbdLayer3D.rotateAllVerts(ones(mesh.N,1), ri, R', COMinitial-linearDx, mesh.p);

    mesh.v = BDF.BDF1(mesh.p, mesh.prevp, h);  
end