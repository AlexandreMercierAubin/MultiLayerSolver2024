% GPU code generation for getting started example (mandelbrot_count.m)
%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.gpuConfig('mex');
cfg.GenerateReport = true;
% cfg.Benchmarking = true;
cfg.GpuConfig.CompilerFlags = '--fmad=false';

maxN = 1e6;
cfg.EnableDynamicMemoryAllocation = true;
cfg.DynamicMemoryAllocationThreshold = 100*maxN;


variableDimension = [false,false,true];
dummyInt = int32(0);
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0.0,[3,3,maxN],variableDimension,'Gpu', true);
ARGS{1}{2} = coder.typeof(0.0,[12,1,maxN],variableDimension,'Gpu', true);

%% Invoke GPU Coder.
codegen -config cfg StVKconstraintGrad3D -args ARGS{1}