function [w] = rigidGeneralizedInverseMass(minv,r,Inertia,n)
%RIGIDGENERALIZEDINVERSEMASS 
%TODO: is minv supposed to be the rigid's entire mass or the vert mass in
%the rigid? I think it is the entire mass.
% r: position wrt center of mass 1x3xnBodies
% Intertia: Intertia matrix 3x3xnBodies
% n : normal of the constraint direction (should be the right line of
% Hstack) 1x3xnBodies

    axis = cross(r,n,2);
    inertiaComponent = pagemldivide(Inertia, pagetranspose(axis));
    w = minv + pagemtimes(axis,inertiaComponent);
end

