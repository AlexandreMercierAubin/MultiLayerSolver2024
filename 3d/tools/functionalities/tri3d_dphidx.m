function dphidx = tri3d_dphidx(V, T)
    assert(size(V,2) == 3);
    assert(size(V,2) == 3);
    numt = size(T,1);
    dphidx = zeros(numt,9);
    for t = 1:numt
        v1Id = T(t,1);
        v2Id = T(t,2);
        v3Id = T(t,3);

        v1 = V(v1Id,:);
        v2 = V(v2Id,:);
        v3 = V(v3Id,:);
        dphidX1 = v2-v1;
        dphidX2 = v3-v1;
        Tlocal = [dphidX1',dphidX2'];
    
        dphidxbot = (Tlocal'*Tlocal)\(Tlocal');
        dphidxlocal(1,:) = -sum(dphidxbot,1);
        dphidxlocal(2:3,:) = dphidxbot;
        dphidx(t,:) = reshape(dphidxlocal,1,9);
    end
end

