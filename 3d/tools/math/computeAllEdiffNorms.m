function [EDiffNorms] = computeAllEdiffNorms(cache,h)
    oldF = cache.oldF;
    if numel(oldF) == 0
        EDiffNorms = zeros(numel(cache.F)/9,1);
        return
    end
    F = cache.F;
    EDiffNorms = mexEdiffNorm3D(F,oldF,h);
end

