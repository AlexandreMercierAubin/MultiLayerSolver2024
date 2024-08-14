function  setRigidTransform( mesh2d, degrees, translation, forceCenterOfMassRotation)
%RIGIDTRANSFORM Transforms mesh points by given rotation and translation.
%Must be called before setRigidMotions
%   mesh.rigidTransform( degrees, translation ) 
    if nargin < 4
        forceCenterOfMassRotation = false;
    end
    c = cosd( degrees );
    s = sind( degrees );
    R = [ c -s; s c ];
    if ~forceCenterOfMassRotation
        p = [ mesh2d.p(1:2:end), mesh2d.p(2:2:end) ]';
        mesh2d.p(1:1:end) = R * p;
    else
        masses = mesh2d.mass(1:2:end);
        mass = sum(masses);
        centerOfMass = [dot(mesh2d.p(1:2:end), masses); dot(mesh2d.p(2:2:end), masses)] / mass;
        p = reshape( mesh2d.p', 2, [] )' - centerOfMass';
        newp = R * p';
        newp = newp' + centerOfMass';
        mesh2d.p = reshape(newp',1,[])';
    end
    
    mesh2d.p(1:2:end) = mesh2d.p(1:2:end) + translation(1);
    mesh2d.p(2:2:end) = mesh2d.p(2:2:end) + translation(2);
    mesh2d.prevp = mesh2d.p;
    mesh2d.prevPrevp = mesh2d.p;
    mesh2d.prevPrevPrevp = mesh2d.p;
end

