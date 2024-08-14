function [B, gamma] = getB(obj, cache)
    % COMPUTEADAPTB Computes B for the rigid-elastic combo.
    %   Here, gamma maps to the FULL set of elastic dofs, but we note that
    %   we only actualy use a subset of the dofs to compute B, i.e., those
    %   DOFs that are adjacent to an elastic element, i.e., elements in the
    %   ActiveBRows list.
    bodyPositions = [obj.RigidBodies.Position]';
    [B, gamma] = obj.peekB(cache,bodyPositions, obj.p);
end

