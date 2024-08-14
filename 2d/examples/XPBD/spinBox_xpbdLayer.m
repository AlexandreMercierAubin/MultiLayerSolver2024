close all;
clear

gravity = 0; % gravity
h = 1/60; % time step

rho = 5;
nu = 0.33; % Poisson ratio: close to 0.5 for rubber
E = 6e4; % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('cantileverP05',resetMesh, material, scale, rot, settings);
mesh2d.setRigidMotion(pi/200,[0,0],h);

mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );
mesh2da3 = AdaptiveMesh( mesh2d );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 10;%using 40 for fair comparison
energyType = 2;
% computeResidual = true;
computeResidual = true;
useMassSplitting = false;
residualType = 1; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = false;

settings.plotFrameTimeHistogram = false;
if settings.plotFrameTimeHistogram
    runUntil = 0.1;%-1 to deactivate
    giveUpEnabled = false;
else
    runUntil = -1;%-1 to deactivate
    giveUpEnabled = false;
end
boundaryCompliance = 1e-9;

integrator = xpbdLayer2D();
% integrator = BackwardEuler();
% integrator = xpbdLayerKinematic2D();
% integrator = xpbd2D();
integrator.Gravity = gravity;
integrator.iterations = [iterations];
% integrator.boundaryCompliance = 0;
% integrator.elasticCompliance = 1e-8;
integrator.computeResiduals = computeResidual;
integrator.energyType = energyType;
integrator.useMassSplitting = useMassSplitting;
integrator.residualType = residualType;
integrator.runUntilResSmallerThan = runUntil;
integrator.useGravityConstraints = useGravityConstraints;

% integrator2 = xpbd2D();
integrator2 = xpbdLayer2D();
integrator2.energyType = energyType;
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = boundaryCompliance;
layers = [100,0];
integrator2.setLayers(layers);
rigidIT = [1];
elasticIterations = iterations - sum(rigidIT);
assert(elasticIterations >= 0);
integrator2.iterations = [rigidIT,elasticIterations];
integrator2.computeResiduals = computeResidual;
integrator2.useMassSplitting = useMassSplitting;
integrator2.giveUpEnabled = giveUpEnabled;
integrator2.LayerOrderingType = layerType;
integrator2.residualType = residualType;
integrator2.runUntilResSmallerThan = runUntil;
integrator2.useGravityConstraints = useGravityConstraints;

integrator3 = xpbdResLayer2D();
integrator3.energyType = energyType;
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = boundaryCompliance;
integrator3.setLayers(layers);
integrator3.iterations = [rigidIT,elasticIterations];
integrator3.computeResiduals = computeResidual;
integrator3.useMassSplitting = useMassSplitting;
integrator3.giveUpEnabled = giveUpEnabled;
integrator3.residualType = residualType;
integrator3.LayerOrderingType = layerType;
integrator3.runUntilResSmallerThan = runUntil;
integrator3.useGravityConstraints = useGravityConstraints;

settings.MakeVideo = 1;
settings.SceneName = 'spinBox_layer';
settings.CamPadding = [1,1,3,2];
settings.campos = [0,1];

settings.RigidificationEnabled = false;
settings.plotLayers = true;

if computeResidual && ~ settings.plotFrameTimeHistogram
    settings.plotResidual = 1;%6, 9, 10, 11
end
settings.residualRange = [-2,-1];
settings.overwriteComparisonPositionFrom1 = true;
% settings.DrawEDots = true;
settings.skipPreProjectionResidual = true;

% td = simulate( {mesh2da,mesh2da2,mesh2da3}, {integrator,integrator2,integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
td = simulate( {mesh2da}, {integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
