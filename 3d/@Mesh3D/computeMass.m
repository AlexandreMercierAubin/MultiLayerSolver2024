function m = computeMass(mesh)
    % computeMass computes the diagonal of the mass matrix
    %   computeMass( mesh, rho ) returns the lumped node diagonal of the mass 
    %   matrix for the given mesh
    if mesh.useAverageMass
        m = zeros(mesh.N * 3, 1);
        rho = cat(1,mesh.materials(mesh.materialIndex).rho);
        for i = 1:size(mesh.t, 1)
            if mesh.elementType(i) == mesh.elementTypeEnum.Shell
                T = mesh.t(i,1:3);
            else
                T = mesh.t(i,:);
            end
            dofInds =[T*3-2;T*3-1;T*3];
            m(dofInds) = m(dofInds) +  rho(i);
        end
        
        for i = 1:mesh.N
            m(i*3-2) = m(i*3-2)/mesh.valence(i);
            m(i*3-1) = m(i*3-2);
            m(i*3) = m(i*3-2);
        end
    else
        %then mass density
        m = zeros(mesh.N * 3, 1);
        rho = cat(1,mesh.materials(mesh.materialIndex).rho);
        for i = 1:size(mesh.t, 1)
            if mesh.elementType(i) == mesh.elementTypeEnum.Shell
                T = mesh.t(i,1:3);
            else
                T = mesh.t(i,:);
            end
            dofInds =[T*3-2;T*3-1;T*3];
            nverts = numel(T);
            % number of verts in the element (tri3D or tet)
            %mass is the density multiplied by the area divided by the
            %number of nodes of the element
            m(dofInds) = m(dofInds) +  (mesh.elV(i) * rho(i))/nverts;
        end
    end
end