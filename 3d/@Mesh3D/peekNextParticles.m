function [p, angularVelocity, rigidBodyPosition, inertia] = peekNextParticles(mesh, h, deltav)
    v = mesh.v;
    v(mesh.activeDOFs) = mesh.v(mesh.activeDOFs) + deltav;
    v(mesh.pinnedDOFs) = v(mesh.pinnedDOFs) * 0;
    p = mesh.p + h * v;
    rigidBodyPosition = [];
    inertia = [];
    angularVelocity = [];
end