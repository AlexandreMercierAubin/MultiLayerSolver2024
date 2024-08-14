function auto_octopus_xpbdLayer(path,percentageImprovement,runOnGPU)
    close all;
    h = 0.01; % time step
    
    % CONFIG
    settings = Simulation3DSettings();
    
    % MESHES
    
    clear meshes;
    %for octopusP2 resolution
    rho = 5;
    nu = 0.45;      % Poisson ratio: close to 0.5 for rubber
    k = 1e7;     % Young's modulus: 0.01e9 approximate for rubber
    [ mu, lambda ] = toLame( nu, k );
    alpha0 = 0.001;   % Rayleigh factor on M
    alpha1 = 0.001;  % Rayleigh factor on K
    tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];
    
    resetMeshCache = true;
    scale = [1,1,1];
    baseMesh = meshLoader("octopus3Fat", [], tMaterial, scale);
    baseMesh.setRigidTransform([90,0,0],[0,0,0],true);
    
    bottom = min(baseMesh.p(3:3:end));
    closeToBottom = baseMesh.p(3:3:end) < bottom + 0.01;
    radius = 0.15;
    nearCenter = abs(baseMesh.p(2:3:end)) < radius & abs(baseMesh.p(1:3:end)) < radius;
    
    pinMiddle = find(closeToBottom & nearCenter);
    baseMesh.pin(pinMiddle);
    
    meshes = AdaptiveMesh3D(baseMesh);
    meshes2 = AdaptiveMesh3D(baseMesh);
    
    gravity = -9.8;
    airDrag = 0.000;
    iterations = 100;
    
    layers = [80,60,40,0];
    rigidIT = [1,2,4];
    iterationsList = [rigidIT,iterations - sum(rigidIT)];
    
    % layers = [0];
    % iterationsList = [iterations];
    
    computeResiduals = true;
    boundaryCompliance = 0;
    
    energyType = 2;
    computeResidual = true;
    residualType = 1; % 0 : constraint residual 1: eq4 residual
    layerType = 0; % 0: strain rate 1: random 2:eigs
    useGravityConstraints = true;
    
    
    settings.plotFrameTimeHistogram = true;
    runUntil = -1;%-1 to deactivate
    giveUpEnabled = true;
    hangingStop = 100;
    settings.StopAtImprovementPercent = percentageImprovement;
    
    
    
    rigidificator = NoRigidificator3D();
    
    integrator = xpbdLayer3D();
    integrator.Gravity = gravity;
    integrator.boundaryCompliance = boundaryCompliance;
    integrator.setLayers(layers);
    integrator.iterations = iterationsList;
    integrator.computeResiduals = computeResidual;
    integrator.giveUpEnabled = giveUpEnabled;
    integrator.residualType = residualType;
    integrator.LayerOrderingType = layerType;
    integrator.runUntilResSmallerThan = runUntil;
    integrator.useGravityConstraints = useGravityConstraints;
    integrator.hangingStop = hangingStop;
    integrator.runOnGPU = runOnGPU;
    integrator.airDrag = airDrag;
    
    integratorElastic = xpbdLayer3D();%elastic
    integratorElastic.Gravity = gravity;
    integratorElastic.boundaryCompliance = boundaryCompliance;
    integratorElastic.contactResistance = 0;
    integratorElastic.setLayers([0]);
    integratorElastic.iterations = [iterations];
    integratorElastic.computeResiduals = computeResidual;
    integratorElastic.giveUpEnabled = giveUpEnabled;
    integratorElastic.residualType = residualType;
    integratorElastic.LayerOrderingType = layerType;
    integratorElastic.runUntilResSmallerThan = runUntil;
    integratorElastic.useGravityConstraints = useGravityConstraints;
    integratorElastic.hangingStop = hangingStop;
    integratorElastic.runOnGPU = runOnGPU;
    integratorElastic.airDrag = airDrag;
    
    energyModel = StVenantKirchoff3DEnergy();
    settings.SceneName = 'auto_octopus_xpbd';
    settings.campos=[-3,-3,1];
    settings.camtarget = [0,0,0];
    settings.residualComparisonLabels = ["E","L"];

    settings.MakeVideo = 1;
    settings.VideoOut = path+string(settings.SceneName)+"/"+string(percentageImprovement);
    settings.FramesToRecord = 6000;
    settings.PlotSkip = 10;
    
    td = simulate3D({meshes, meshes2},h,NullContactFinder(), {integratorElastic,integrator}, rigidificator, settings, energyModel,NullAnimationScript());
    save(path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName)+".mat", 'td');
    writeTDcsv(td, path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName), ["_elastic", "_xpbdResLayer"]);
end