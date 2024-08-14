cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

rho = 1;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 5e2;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu,k );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.01;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

baseMesh = polyLoader("cube",  tMaterial,[1,1,1], true, settings);

baseMesh.setRigidTransform([90,0,0],[0,0,0.4]);

meshes = AdaptiveMesh3D(baseMesh);

gravity = -9.8;
layers = [50,25,0];
rigidIT = [1,4];

iterations = 20;

contactCompliance = 10^-4;
computeResiduals = true;
contactResistance = 0.90;
boundaryCompliance = 0;

energyType = 2;
computeResidual = true;
useMassSplitting = false;
residualType = 1; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;

hangingStop = -1;
settings.plotFrameTimeHistogram = false;
if settings.plotFrameTimeHistogram
    runUntil = 0.5;%-1 to deactivate
    giveUpEnabled = false;
    hangingStop = 500;
else
    runUntil = -1;%-1 to deactivate
    giveUpEnabled = false;
end

rigidificator = NoRigidificator3D();

integrator = xpbdLayer3D();
integrator.Gravity = gravity;
integrator.boundaryCompliance = boundaryCompliance;
integrator.contactResistance = contactResistance;
integrator.setLayers(layers);
integrator.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator.computeResiduals = computeResidual;
integrator.useMassSplitting = useMassSplitting;
integrator.giveUpEnabled = giveUpEnabled;
integrator.residualType = residualType;
integrator.LayerOrderingType = layerType;
integrator.runUntilResSmallerThan = runUntil;
integrator.useGravityConstraints = useGravityConstraints;
integrator.hangingStop = hangingStop;

energyModel = NeoHookean3DEnergy();

settings.MakeVideo = 1;
settings.SceneName = 'cubeDrop';
% settings.FramesToRecord = 600;
% settings.DrawRigidDOFs = true;
settings.PlotEDotHist = 1;
settings.campos=[7,7,1];
% settings.quicksolveSimulation = 1;
% settings.PCGiterations = 100;
% settings.PGSiterations = 100;
settings.recomputeCacheAinv = true;

simulate3D(meshes,h,NullContactFinder(), integrator, rigidificator, settings, energyModel);