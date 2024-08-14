function [activeDOFs, ActiveElasticDOFs, ActiveDofsCorrespondingID] = computeActiveDOFs(ElasticDOFs, unpinnedDOFs, totalDOFs, isRigidComponentPinned)
    %COMPUTEACTIVEDOFS Finds dofs that are both unpinned and elastic
    dofIdRange = 1:totalDOFs;
    setDiffIds = true(numel(dofIdRange),1);
    setDiffIds(ElasticDOFs) = false;
    rigid = dofIdRange(setDiffIds);

    logical_idx = false(1,totalDOFs);
    logical_idx(unpinnedDOFs) = true;
    scaledLogical = logical_idx;
    scaledLogical(rigid) = false;
    logical_idx(rigid) = [];
    
    ActiveDofsCorrespondingID = scaledLogical;
    ActiveElasticDOFs = find(logical_idx);
    elasticN = numel(ElasticDOFs);
    if numel(isRigidComponentPinned) < 1
        activeDOFs = ActiveElasticDOFs;
        return;
    end
    isRigidUnpinned = find(~isRigidComponentPinned);
    rigidUnpinned = elasticN+reshape([isRigidUnpinned*3-2;isRigidUnpinned*3-1;isRigidUnpinned*3],1,[]);

    activeDOFs = [ActiveElasticDOFs, rigidUnpinned];
end