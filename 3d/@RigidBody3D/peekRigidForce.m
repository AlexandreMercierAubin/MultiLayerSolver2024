function [Force, Torque] = peekRigidForce(body, h, AngularVelocity, RigidBodyPosition, Inertia, meshf, meshp)
    %COMPUTERIGIDFORCE Computes the force and torque on a rigid body using
    %the force of the individual particles
    
    Force = body.Force + sum(reshape(meshf(body.DOFs), 3, []), 2);
    
    p = reshape(meshp(body.DOFs), 3, []);
    ri = p - RigidBodyPosition;
    f = reshape(meshf(body.DOFs), 3, []);
    torque = cross(ri, f);
    
    torqueSum = sum(torque,2);
    Torque = body.Torque + torqueSum;
    
    % And also the omega cross J omega term...
    
    L = Inertia*AngularVelocity;
    Lhat = crossProductMatrix(L)*h;

    Jmod = Inertia - h*h*Lhat*diag(1./[body.Mass,body.Mass,body.Mass])*Lhat;
    
    Torque = Torque - cross( AngularVelocity, Jmod * AngularVelocity );
end
