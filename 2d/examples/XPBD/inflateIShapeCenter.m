close all;
clear

gravity = 0; % gravity
h = 0.0025; % time step

rho = 1;
nu = 0.36; % Poisson ratio: close to 0.5 for rubber
E = 1e3; % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.0001; % Rayleigh factor on M
alpha1 = 0.001; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('IShapeHD.1',resetMesh, material, scale, rot, settings);

% inds = find(mesh2d.p(2:2:end)< 0.7 & mesh2d.p(2:2:end) > 0.5);
% mesh2d.p(inds*2) = 1.07*mesh2d.p(inds*2);
mesh2d.p = 1.12*mesh2d.p;

mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );
mesh2da3 = AdaptiveMesh( mesh2d );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 50;
energyType = 2;
computeResidual = true;
useMassSplitting = false;
residualType = 1; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;
if residualType == 1
    useGravityConstraints = false;
end

settings.plotFrameTimeHistogram = true;
hangingStop = -1;
if settings.plotFrameTimeHistogram
    runUntil = 0.2;%-1 to deactivate
    giveUpEnabled = false;
    hangingStop = 80;
else
    runUntil = -1;%-1 to deactivate
    giveUpEnabled = false;
end
% settings.residualHistBinSize = 0.0025;
boundaryCompliance = 1e-9;

integrator = xpbdLayer2D();
integrator.Gravity = gravity;
integrator.iterations = [iterations];
integrator.computeResiduals = computeResidual;
integrator.energyType = energyType;
integrator.useMassSplitting = useMassSplitting;
integrator.residualType = residualType;
integrator.runUntilResSmallerThan = runUntil;
integrator.useGravityConstraints = useGravityConstraints;
integrator.hangingStop = hangingStop;

integrator2 = xpbdLayer2D();
integrator2.energyType = energyType;
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = boundaryCompliance;
layers = [80,60,0];
integrator2.setLayers(layers);
rigidIT = [5,3];
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
integrator2.hangingStop = hangingStop;

integrator3 = xpbdResLayer2D();
integrator3.energyType = energyType;
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = boundaryCompliance;
integrator3.setLayers(layers);
integrator3.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator3.computeResiduals = computeResidual;
integrator3.useMassSplitting = useMassSplitting;
integrator3.giveUpEnabled = giveUpEnabled;
integrator3.residualType = residualType;
integrator3.LayerOrderingType = layerType;
integrator3.runUntilResSmallerThan = runUntil;
integrator3.useGravityConstraints = useGravityConstraints;
integrator3.hangingStop = hangingStop;

settings.MakeVideo = 1;
settings.SceneName = 'inflateIShape_layer';
settings.CamPadding(4) = 1;
settings.CamPadding(3) = 2;
settings.campos = [0,1];

settings.RigidificationEnabled = false;
settings.plotLayers = true;
settings.residualRange = [-1,0.2];
% settings.DrawEDots = true;
settings.overwriteComparisonPositionFrom1 = true;
settings.FramesToRecord = 0.15/h;
settings.InitialWindowPosition = [0, 0, 1920, 1080];
if computeResidual && ~settings.plotFrameTimeHistogram
    settings.plotResidual = 1;%6, 9, 10, 11
end

td = simulate( {mesh2da,mesh2da2,mesh2da3}, {integrator,integrator2,integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2da,mesh2da2}, {integrator,integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
writeTDcsv(td, string(settings.SceneName), ["_elastic","_xpbdLayer", "_xpbdResLayer"]);
% td = simulate( {mesh2da}, {integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
