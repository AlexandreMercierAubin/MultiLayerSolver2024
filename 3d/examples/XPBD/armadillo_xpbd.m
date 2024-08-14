cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear mesha;

rho = 5;
nu = 0.37;      % Poisson ratio: close to 0.5 for rubber
k = 1e6;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu,k );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
tMaterialArmadillo = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];
tMaterialArmadillo.cacheName = 'armadilloSingleMat';

settings.MakeVideo = 1;
settings.SceneName = 'armadillo_xpbd';
% settings.FramesToRecord = 500;
% settings.InitialWindowPosition = [0,0,1920,1080];
settings.campos=[0,10,3];

resetMesh = true;
mesh = polyLoader("arma_6", tMaterialArmadillo,[1,1,1],resetMesh, settings);

mesh.setRigidTransform([90,0,0],[0,0,0.1]);
mesha = AdaptiveMesh3D(mesh);

gravity = -9.8;
% gravity = 0;
layers = [100,50,30,0];
rigidIT = [1,2,2];
% layers = [100,0];
% rigidIT = [1];

iterations = 30;

computeResiduals = true;
contactResistance = 0.5;
boundaryCompliance = 0;

energyType = 2;
computeResidual = true;
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

integrator3 = xpbdLayer3D();
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = boundaryCompliance;
integrator3.contactResistance = contactResistance;
integrator3.setLayers(layers);
integrator3.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator3.computeResiduals = computeResidual;
integrator3.giveUpEnabled = giveUpEnabled;
integrator3.residualType = residualType;
integrator3.LayerOrderingType = layerType;
integrator3.runUntilResSmallerThan = runUntil;
integrator3.useGravityConstraints = useGravityConstraints;
integrator3.hangingStop = hangingStop;

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,0], 0);
meshMeshContactFinder = MeshSCD3D(0.0,false, []);
contactFinder = {planeContactFinder, meshMeshContactFinder};

energyModel = NeoHookean3DEnergy();

simulate3D(mesha,h,contactFinder, integrator3, rigidificator, settings, energyModel);