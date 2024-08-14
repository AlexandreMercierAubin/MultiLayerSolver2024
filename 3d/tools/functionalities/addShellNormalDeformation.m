function addShellNormalDeformation(meshes, p, cache,settings)
    if settings.addShellNormalDeformation
        % we need to add the normal change to the deformation gradient of shells
        numt = size(meshes.t,1);
        n_vec = zeros(9*numt,1);
        cache.dndp = zeros(3,9,numt);
        localN = zeros(9,3);

        %preallocates the sparse matrix with the biggest case scenario in
        %mind (all shells)... can be improved with the number of shell
        %elements instead
        ii = zeros(1,9*9*numt);
        jj = zeros(1,9*9*numt);
        vals = zeros(1,9*9*numt);
        counter = 0;

        for i = 1:numt
            t = i;
            if meshes.elementType(t) ~= meshes.elementTypeEnum.Shell
                continue;
            end

            v1Id = meshes.t(t,1);
            v1 = p(v1Id*3-2:v1Id*3);
            v2Id = meshes.t(t,2);
            v2 = p(v2Id*3-2:v2Id*3);
            v3Id = meshes.t(t,3);
            v3 = p(v3Id*3-2:v3Id*3);
            
            localN(1:3,1) = meshes.referenceSpaceNormals(i,:);
            localN(4:6,2) = meshes.referenceSpaceNormals(i,:);
            localN(7:9,3) = meshes.referenceSpaceNormals(i,:);

            x = [v1;v2;v3];
            I1 = [-eye(3), eye(3), zeros(3)];
            I2 = [-eye(3), zeros(3), eye(3)];
            laplacex1 = I1 * x;
            laplacex2 = I2 * x;
            crossx1 = crossProductMatrix(laplacex1);
            crossx2 = crossProductMatrix(laplacex2);
            ntilde = cross(laplacex1,laplacex2);
            normntilde = norm(ntilde);
            assert(~isinf(normntilde));%if this triggers, then the normal is made of zeros... check your mesh
            localn = ntilde./normntilde;          
            localnvec =  localN*localn;
            n_vec(t*9-8:t*9) = localnvec;
            
            invNtilde = 1/normntilde;

            dndp = invNtilde*(eye(3)-localn*localn')*(crossx1*I2 - crossx2*I1); 
            cache.dndp(:,:,t) = dndp;
            Ndndp = localN*dndp;

            counter = counter + 1;
            idi = repmat(t*9-8:t*9,1,9);
            ii(counter*81-80:counter*81) = idi;
            idj1 = reshape(repmat(meshes.t(t,1)*3-2:meshes.t(t,1)*3,9,1),[],1)';
            idj2 = reshape(repmat(meshes.t(t,2)*3-2:meshes.t(t,2)*3,9,1),[],1)';
            idj3 = reshape(repmat(meshes.t(t,3)*3-2:meshes.t(t,3)*3,9,1),[],1)';
            jj(counter*81-80:counter*81) = [idj1,idj2,idj3];
            B1 = reshape(Ndndp(1:9,1:3),1,[]);
            B2 = reshape(Ndndp(1:9,4:6),1,[]);
            B3 = reshape(Ndndp(1:9,7:9),1,[]);
            vals(counter*81-80:counter*81) = [B1,B2,B3];
        end
        %scales back the vector to the actual number of entries
        ii = ii(1:counter*81);
        jj = jj(1:counter*81);
        
        vals = vals(1:counter*81);

        fullNdndp = sparse(ii,jj,vals,9*size(meshes.t,1),numel(p));
        meshes.B = meshes.BnoNormal + fullNdndp;
        cache.F = meshes.B*p + n_vec;
    else
        cache.F = meshes.B*p;
    end
end