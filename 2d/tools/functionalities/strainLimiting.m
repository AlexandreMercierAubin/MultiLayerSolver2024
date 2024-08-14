function strainLimiting(mesh2D,cache)
% strainLimiting(mesh2D,cache)
% applies strain limits from the materials to the mesh
%based on http://graphics.berkeley.edu/papers/Wang-MRI-2010-12/Wang-MRI-2010-12.pdf
    numt = size(mesh2D.t,1);
    for i = 1:numt
        if isa(mesh2D,"AdaptiveMesh") && ~mesh2D.isTriElastic(i)
            continue;
        end
        localF = reshape(cache.F(i*4-3:i*4),2,2);
        [U,S,V] = svd(localF);
        matId = mesh2D.materialIndex(i);
        upperbound = mesh2D.materials(matId).strainUpperBound;
        lowerbound = mesh2D.materials(matId).strainLowerBound;
        Sdiag = diag(S);
        Sdiag(Sdiag > upperbound) = upperbound;
        Sdiag(Sdiag < lowerbound) = lowerbound;

        localFStar = U*diag(Sdiag)*V';
        %strain limiting inversion safeguard
        count = 0;
        while det(localFStar) < 0 && count <2
            minSdiag = min(Sdiag);
            Sdiag(Sdiag == minSdiag) = -minSdiag;
            localFStar = U*diag(Sdiag)*V';
            count = count + 1;
        end

        cache.F(i*4-3:i*4) = reshape(localFStar,4,1);

        Dx = localFStar*mesh2D.el(i).restLength;

        triangleVert = mesh2D.t(i,:);
        ids = [triangleVert*2-1;triangleVert*2];

        %center of mass
        mass = mesh2D.mass(triangleVert*2);
        sumMass = sum(mass);
        d1x = Dx(:,1);
        d2x = Dx(:,2);
        currentPos = mesh2D.p(ids);
        massWeightedPos = currentPos*mass;
        cx = massWeightedPos ./ sumMass;
        DxMass = (d1x*mass(2)+d2x*mass(3));
        WeightedDx = DxMass./sumMass;
        x0 = cx - WeightedDx;
        x1 = x0 + d1x;
        x2 = x0 + d2x;
        mesh2D.p(triangleVert(1)*2-1:triangleVert(1)*2)=x0;
        mesh2D.p(triangleVert(2)*2-1:triangleVert(2)*2)=x1;
        mesh2D.p(triangleVert(3)*2-1:triangleVert(3)*2)=x2;
    end
end

