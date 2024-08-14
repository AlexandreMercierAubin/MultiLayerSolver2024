function [Force,Torque] = peekRigidForce(body,f)
    %COMPUTERIGIDFORCE Computes the force and torque on a rigid body using
    %the force of the individual particles
    Force = body.Force + sum(reshape(f(body.DOFs), 2, []), 2);

    pxfy = (body.Mesh.p(body.Indices * 2 - 1) - body.Position(1)) .* f(body.Indices * 2);
    pyfx = (body.Mesh.p(body.Indices * 2)     - body.Position(2)) .* f(body.Indices * 2 - 1);
    
    Torque = sum(pxfy - pyfx); % fast cross product
end
