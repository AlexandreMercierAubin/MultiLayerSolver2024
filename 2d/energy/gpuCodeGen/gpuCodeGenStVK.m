% GPU code generation for getting started example (mandelbrot_count.m)
%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.gpuConfig('mex');
cfg.GenerateReport = true;
% cfg.Benchmarking = true;
cfg.GpuConfig.CompilerFlags = '--fmad=false';

dummyInt = int32(0);
ARGS = cell(1,1);
ARGS{1} = cell(9,1);
ARGS{1}{1} = coder.typeof(0,[inf,1]);
ARGS{1}{2} = coder.typeof(0,[inf,3]);
ARGS{1}{3} = coder.typeof(dummyInt,[inf,6]);
ARGS{1}{4} = coder.typeof(0,[inf,1]);
ARGS{1}{5} = coder.typeof(0,[inf,6]);
ARGS{1}{6} = coder.typeof(0,[3,3,inf]);
ARGS{1}{7} = coder.typeof(0,[inf,1]);
ARGS{1}{8} = coder.typeof(0,[2,2,inf]);
ARGS{1}{9} = coder.typeof(0,[1,1]);

%% Invoke GPU Coder.
codegen -config cfg StVKForColor -args ARGS{1}