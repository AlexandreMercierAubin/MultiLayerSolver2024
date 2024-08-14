function [p, angularVelocity, rigidBodyPosition, inertia] = peekNextParticles(mesh, h, deltav, ActiveElasticDOFs, ActiveDofsCorrespondingID)
    %UPDATEPARTICLES updates the particles of the mesh given a delta v
    %vector
    elasticN = numel(ActiveElasticDOFs);

    v = mesh.v;
    v(mesh.ActiveDofsCorrespondingID) = v(ActiveDofsCorrespondingID) + deltav(1:elasticN);
    p = mesh.p;
    p(mesh.ActiveDofsCorrespondingID) = mesh.p(ActiveDofsCorrespondingID) + h * v(ActiveDofsCorrespondingID); 

    index = 1;

    angularVelocity = zeros(3,numel(mesh.RigidBodies));
    rigidBodyPosition = zeros(3,numel(mesh.RigidBodies));
    inertia = zeros(3,3,numel(mesh.RigidBodies));
    for i = 1:numel(mesh.RigidBodies)
        rigidBodyPosition(:,i) = mesh.RigidBodies(i).Position;
        if mesh.RigidBodies(i).isPinned
            mesh.v(mesh.RigidBodies(i).DOFs) = 0;
            continue;
        end
        [pList,vList,inds, bodyp, av, bodyInertia] = mesh.RigidBodies(i).peekPosition(h, deltav(elasticN + index * 6 - 5:elasticN + index * 6), p);
        rigidBodyPosition(:,i) = bodyp;
        angularVelocity(:,i) = av;
        inertia(:,:,i) = bodyInertia;
        p(inds) = pList;
        v(inds) = vList;
        index = index + 1;
    end
end