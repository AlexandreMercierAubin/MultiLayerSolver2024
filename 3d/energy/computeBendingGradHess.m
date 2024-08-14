function [psi, f, shellBendingH]=computeBendingGradHess(mesh3D, p, cache, settings)
    psi = 0;
    if nargout > 1
        f = zeros(mesh3D.N*3,1);
	    shellBendingH = zeros(mesh3D.N*3,mesh3D.N*3);
    end

	if settings.addBendingEnergy

        %TODO: use the normals from addShellNormalDeformation
        V = mesh3D.formatPositions(p);
        triNormalDirections = normals(V,mesh3D.t,"Stable",true);  %it computes more normals than needed, but whatever
        %gptoolbox does not normalize normals so we do it here
        triNormals = normr(triNormalDirections);

        innerAngles = internalangles(V,mesh3D.t); %a23 a13 a12 -> v1 v2 v3

        areas = triangleArea3D(V,mesh3D.t);

        norm1 = triNormals(mesh3D.bendingEdgesElements(:,1),:);
        norm2 = triNormals(mesh3D.bendingEdgesElements(:,2),:);
        %get the axis of rotation
        vertex1full = mesh3D.edgeVertices(:,1);
        vertex2full = mesh3D.edgeVertices(:,2);
        v1full = p([vertex1full*3-2,vertex1full*3-1,vertex1full*3]);
        v2full = p([vertex2full*3-2,vertex2full*3-1,vertex2full*3]);
        e0full = v2full-v1full;
        e0norm = normr(e0full);
        angles = mexAngleBetweenVectors(norm1,norm2,e0norm);
        anglesDiff = angles-mesh3D.edgeRestDihedral;
        kd = mesh3D.kdRestEdgeAverageHeight;
        psi_i = kd.*anglesDiff.*anglesDiff;
        psi = sum(psi_i);

        %TODO: probably needs to store the triangle edges for bending
        for i = 1:size(mesh3D.bendingEdges,1)
            triIDs = mesh3D.bendingEdgesElements(i,:);
            vertex1 = mesh3D.edgeVertices(i,1);
            vertex2 = mesh3D.edgeVertices(i,2);
            opposite1 = mesh3D.edgeOppositeVertices(i,1);
            opposite2 = mesh3D.edgeOppositeVertices(i,2);
            
            %finds the right order for the innerAngles
            t1 = mesh3D.t(triIDs(1),:);
            t2 = mesh3D.t(triIDs(2),:);
            for j = 1:3
                if t1(j) == vertex1
                    innerAngle1 = innerAngles(triIDs(1),j);
                end
                if t2(j) == vertex1
                    innerAngleTilde1 = innerAngles(triIDs(2),j);
                end

                if t1(j) == vertex2
                    innerAngle2 = innerAngles(triIDs(1),j);
                end
                if t2(j) == vertex2
                    innerAngleTilde2 = innerAngles(triIDs(2),j);
                end
            end

            v1 = p(vertex1*3-2:vertex1*3);
            v2 = p(vertex2*3-2:vertex2*3);
            o1 = p(opposite1*3-2:opposite1*3);
            o2 = p(opposite2*3-2:opposite2*3);

            area = areas(triIDs(1));
            areaTilde = areas(triIDs(2));

            ids = [opposite1*3-2:opposite1*3,vertex1*3-2:vertex1*3,vertex2*3-2:vertex2*3,opposite2*3-2:opposite2*3];
            if nargout > 1
                [g,H] = neoGrinspunianBendingEnergy(kd(i),anglesDiff(i),[o1';v1';v2';o2'],norm1(i,:),norm2(i,:),area,areaTilde,innerAngle1,innerAngleTilde1,innerAngle2,innerAngleTilde2);

                shellBendingH(ids,ids) = shellBendingH(ids,ids) + H;

                f(ids) = f(ids) + g;
            end
        end

        if nargout > 1
            shellBendingH = sparse(shellBendingH); %todo: make this directly into a sparse instead of dense to sparse.
        end
        %-------------------------------------------------------------
        %mexed version

%         [ii,jj,Hvals, forces, angles, psi] = mexGrinspunBendingGradHess(p, mesh3D.bendingEdges,mesh3D.edgeOppositeVertices, mesh3D.kdRestEdgeAverageOppositeHeight, mesh3D.edgeRestDihedral);
%         cache.shellBendingH = sparse(ii,jj,Hvals, mesh3D.N*3, mesh3D.N*3);
%         cache.elasticForces = cache.elasticForces + forces;
    end
end

