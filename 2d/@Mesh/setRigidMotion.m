function setRigidMotion( mesh, angularVelocity, velocity, h, debugplot)
%SETRIGIDMOTION Sets mesh point velocities to have given velocity
%   mesh.rigidTransform( radPerSec, velocity )
%   Angular motion is about the center of mass

    if nargin < 5
        debugplot = false;
    end
    [~, COMinitial, angle, ~, ~, ri, ~, ~]=AdaptiveMesh.makeRigidsFromConnectivity(1, mesh.p, zeros(size(mesh.p)), ones(mesh.N,1), mesh.mass);
    [R, COM, ~] = RigidBody.rigidSymplecticImperativeStep(COMinitial, angle, h, velocity, angularVelocity);
    xi = rotateAllVerts2D(ones(mesh.N,1), ri, R, COM, mesh.p);
    
    linearDx = h*velocity;
    mesh.prevp = rotateAllVerts2D(ones(mesh.N,1),ri,R',COMinitial-linearDx,mesh.p);
    mesh.prevPrevp = rotateAllVerts2D(ones(mesh.N,1),ri,R'*R',COMinitial-2*linearDx,mesh.p);
    mesh.prevPrevPrevp = rotateAllVerts2D(ones(mesh.N,1),ri,R'*R'*R',COMinitial-3*linearDx,mesh.p);

    %it is better to only add the linear velocities after the rotation
    % v = BDF.BDF3(mesh.p, mesh.prevp, mesh.prevPrevp, mesh.prevPrevPrevp, h);
    % v = BDF.BDF2(mesh.p, mesh.prevp, mesh.prevPrevp, h);
    v = BDF.BDF1(mesh.p, mesh.prevp, h);%

    mesh.v = v;

    if debugplot
        hold on;
        debugDot(COM,'red',20);
        debugDot(COMinitial,'yellow',20)
        debugPlotVectors(xi,mesh.v-repmat(velocity',mesh.N,1),"green");
        % debugPlotVectors(xi,mesh.v,"blue");
    end
end