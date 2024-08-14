% GPU code generation for getting started example (mandelbrot_count.m)
%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.gpuConfig('mex');
cfg.GenerateReport = true;
% cfg.Benchmarking = true;
cfg.GpuConfig.CompilerFlags = '--fmad=false';

%first param is just a type so it enforces type safety
dummyInt = int32(0);
dummyLogical = true;
dummyDouble = 0.0;
ARGS = cell(1,1);
ARGS{1} = cell(44,1);

ARGS{1}{1} = coder.typeof(dummyLogical,[inf,inf]); %isElasticElement
ARGS{1}{2} = coder.typeof(dummyLogical,[inf,inf]); %isRigidElement
ARGS{1}{3} = coder.typeof(dummyInt,[inf,inf]); %rigidBodySetsPerLayer
ARGS{1}{4} = coder.typeof(dummyInt,[inf,inf]); %vertexRigidBody
ARGS{1}{5} = coder.typeof(dummyLogical,[inf,inf]); %isBoundaryVertex
ARGS{1}{6} = coder.typeof(dummyLogical,[inf,inf]); %isRigidVertex
ARGS{1}{7} = coder.typeof(dummyDouble,[inf,1]); %xi
ARGS{1}{8} = coder.typeof(dummyDouble,[inf,3]); %lambdai
ARGS{1}{9} = coder.typeof(dummyDouble,[1,1]); %sumNonMonitoredTime
ARGS{1}{10} = coder.typeof(dummyDouble,[1,1]); %h
ARGS{1}{11} = coder.typeof(dummyDouble,[inf,1]); %numRigids
ARGS{1}{12} = coder.typeof(dummyDouble,[inf,inf]); %residualArray
ARGS{1}{13} = coder.typeof(dummyDouble,[1,1]); %boundaryCompliance
ARGS{1}{14} = coder.typeof(dummyDouble,[1,1]); %contactResistance
ARGS{1}{15} = coder.typeof(dummyDouble,[1,1]); %contactCompliance
ARGS{1}{16} = coder.typeof(dummyDouble,[1,inf]); %iterations
ARGS{1}{17} = coder.typeof(dummyDouble,[1,inf]); %layers
ARGS{1}{18} = coder.typeof(dummyDouble,[1,1]); %runUntilResSmallerThan
ARGS{1}{19} = coder.typeof(dummyLogical,[1,1]); %computeResiduals
ARGS{1}{20} = coder.typeof(dummyDouble,[1,1]); %hangingStop
ARGS{1}{21} = coder.typeof(dummyLogical,[1,1]); %giveUpEnabled
ARGS{1}{22} = coder.typeof(dummyDouble,[1,1]); %giveUpThreshold
ARGS{1}{23} = coder.typeof(dummyLogical,[1,1]); %LastLayerGiveUpEnabled
ARGS{1}{24} = coder.typeof(dummyLogical,[1,1]); %useGravityConstraints
ARGS{1}{25} = coder.typeof(dummyDouble,[inf,1]); %oldp
ARGS{1}{26} = coder.typeof(dummyDouble,[inf,1]); %oldOldp
ARGS{1}{27} = coder.typeof(dummyDouble,[inf,2]); %contactNormals
ARGS{1}{28} = coder.typeof(dummyDouble,[inf,2]); %pointColliders
ARGS{1}{29} = coder.typeof(dummyDouble,[inf,1]); %contactVertexID
ARGS{1}{30} = coder.typeof(dummyDouble,[inf,1]); %mass
ARGS{1}{31} = coder.typeof(dummyInt,[inf,3]); %T
ARGS{1}{32} = coder.typeof(dummyDouble,[1,inf]); %unpinnedDOFs
ARGS{1}{33} = coder.typeof(dummyDouble,[1,inf]); %pinnedDOFs
ARGS{1}{34} = coder.typeof(dummyDouble,[2,2,inf]); %DmInv
ARGS{1}{35} = coder.typeof(dummyDouble,[3,3,inf]); %perElementXPBDalphaMatrix
ARGS{1}{36} = coder.typeof(dummyDouble,[3,3,inf]); %perElementXPBDalphaInvMatrix
ARGS{1}{37} = coder.typeof(dummyDouble,[inf,1]); %elementColor
ARGS{1}{38} = coder.typeof(dummyInt,[inf,6]); %elementDOFs
ARGS{1}{39} = coder.typeof(dummyDouble,[1,1]); %numColors
ARGS{1}{40} = coder.typeof(dummyDouble,[inf,1]); %f
ARGS{1}{41} = coder.typeof(dummyDouble,[inf,1]); %layermap
ARGS{1}{42} = coder.typeof(dummyDouble,[1,1]); %gravityCompliance
ARGS{1}{43} = coder.typeof(dummyDouble,[1,1]); %gravity
ARGS{1}{44} = coder.typeof(dummyDouble,[3,3,inf]);%perElementXPBDBeta

%% Invoke GPU Coder.
codegen -config cfg wrapLayerSolve -args ARGS{1}