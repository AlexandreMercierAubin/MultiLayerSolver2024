function [EDotNorms] = computeAllEdotNorms(mesh3D, cache)
    F = mesh3D.B * (mesh3D.p + cache.oldp)/2;
    Fdot = mesh3D.B * mesh3D.v;
    
    EDotNorms = mexEdotNorm3D( F, Fdot );
end

