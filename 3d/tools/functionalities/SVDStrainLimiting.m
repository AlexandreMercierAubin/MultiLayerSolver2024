function SVDStrainLimiting(meshes, p, cache, settings)
    %SVDStrainLimiting(meshes, cache, settings)
    %using material properties, limits strain by SVD decomposition and
    %analyzing eigenvalues
    
    %TODO: this can cause pinned dofs to move due to strain limits. Make
    %sure that pinned dofs don't move
    if settings.StrainLimitingEnabled
        %based on http://graphics.berkeley.edu/papers/Wang-MRI-2010-12/Wang-MRI-2010-12.pdf
        numt = size(meshes.t,1);
        for i = 1:numt
            if isa(meshes,"AdaptiveMesh3D") && ~meshes.isTetElastic(i)
                continue;
            end
            localF = reshape(cache.F(i*9-8:i*9),3,3);
            [U,S,V] = svd(localF);
            matId = meshes.materialIndex(i);
            upperbound = meshes.materials(matId).strainUpperBound;
            lowerbound = meshes.materials(matId).strainLowerBound;

            %TODO:figure out why S is reduced on rigidification.
            %clamp the eigenvalues of F
            Sdiag = diag(S);
            Sdiag(Sdiag > upperbound) = upperbound;
            Sdiag(Sdiag < lowerbound) = lowerbound;


            %strain limiting inversion safeguard
            %should this be before clamping?
            if det(localF) < 0 
                minSdiag = min(Sdiag);
                Sdiag(Sdiag == minSdiag) = -minSdiag;
            end

            %reconstructing a corrected deformation gradient
            localFStar = U*diag(Sdiag)*V';
            cache.F(i*9-8:i*9) = reshape(localFStar,9,1);

            triangleVert = meshes.t(i,:);
            Dr = meshes.restFrame(:,:,i);
            if meshes.elementType == meshes.elementTypeEnum.Shell
                Dr = Dr(:,1:2);
                triangleVert = triangleVert(:,1:3);
            end
            pinnedVerts = meshes.pinned(triangleVert);
            pinnedDofs = logical(reshape(repmat(pinnedVerts,1,3)',[],1));
            ids = [triangleVert*3-2; triangleVert*3-1;triangleVert*3];
            currentPos = p(ids);

            %TODO: figure out why this F needs to be transposed for the
            %rotations to work properly...
            Dx = localFStar'*Dr;

            %center of mass
            mass = meshes.mass(triangleVert*3);%mass is repeated per dofs so 3x N particles
            massEnd = mass(2:end);

            sumMass = sum(mass);
            DxMass = Dx*massEnd;
            weightedDistance = DxMass./sumMass;

            
            weightedPos = currentPos * mass;
            centerOfMass = weightedPos ./ sumMass;
            x0 = centerOfMass - weightedDistance;

            xn = [x0,x0+Dx];
            xn(pinnedDofs) = currentPos(pinnedDofs);
            p(ids) = xn;
        end
    end
end