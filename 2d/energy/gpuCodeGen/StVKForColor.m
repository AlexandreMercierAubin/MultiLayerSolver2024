function [lambdaParallel,xiParallel] = StVKForColor(lambdaParallel, xiParallel, xiPrevParallel, alphaTilde, minvParallel, DmInv, betaTilde, h, volume)
    %generate stvk energy on gpu with cuda
    
    % coder.extrinsic("mexStVKVoigtConstraintEigen");
    % Add Kernelfun pragma to trigger kernel creation
    coder.gpu.kernelfun(); %https://www.mathworks.com/help/gpucoder/ref/coder.gpu.kernel.html
    [deltaxi, ~, deltaLambdai, ~] = matStvkVoigtConstraint( xiParallel, xiPrevParallel, DmInv, alphaTilde, lambdaParallel, minvParallel, betaTilde, h, volume);
    lambdaParallel = lambdaParallel + deltaLambdai;
    xiParallel = xiParallel + deltaxi;
end