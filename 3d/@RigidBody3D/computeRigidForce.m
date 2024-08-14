function computeRigidForce(body, h)
    %COMPUTERIGIDFORCE Computes the force and torque on a rigid body using
    %the force of the individual particles
    [body.Force, body.Torque] = peekRigidForce(body, h, body.AngularVelocity, body.Position, body.Inertia, body.Mesh.f, body.Mesh.p);
end
