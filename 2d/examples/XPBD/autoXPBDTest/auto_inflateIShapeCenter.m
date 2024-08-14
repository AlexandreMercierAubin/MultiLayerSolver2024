function auto_inflateIShapeCenter(path,percentageImprovement,runOnGPU,binSize)
close all;

gravity = 0; % gravity
h = 0.0025; % time step

rho = 3;
nu = 0.35; % Poisson ratio: close to 0.5 for rubber
E = 4e6; % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.0001; % Rayleigh factor on M
alpha1 = 0.001; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();
runUntil = 0;
settings.StopAtImprovementPercent = percentageImprovement;

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

iterations = 100;
energyType = 2;
computeResidual = true;
useMassSplitting = false;
residualType = 0; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs

settings.plotFrameTimeHistogram = true;
giveUpEnabled = false;
hangingStop = 100;
boundaryCompliance = 1e-9;

integrator = xpbdLayer2D();
integrator.Gravity = gravity;
integrator.iterations = [iterations];
integrator.computeResiduals = computeResidual;
integrator.energyType = energyType;
integrator.useMassSplitting = useMassSplitting;
integrator.residualType = residualType;
integrator.runUntilResSmallerThan = runUntil;
integrator.hangingStop = hangingStop;
integrator.runOnGPU = runOnGPU;

integrator2 = xpbdLayer2D();
integrator2.energyType = energyType;
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = boundaryCompliance;
layers = [100,66,33,0];
integrator2.setLayers(layers);
rigidIT = [1,2,4];
elasticIterations = iterations - sum(rigidIT);
assert(elasticIterations >= 0);
integrator2.iterations = [rigidIT,elasticIterations];
integrator2.computeResiduals = computeResidual;
integrator2.useMassSplitting = useMassSplitting;
integrator2.giveUpEnabled = giveUpEnabled;
integrator2.LayerOrderingType = layerType;
integrator2.residualType = residualType;
integrator2.runUntilResSmallerThan = runUntil;
integrator2.hangingStop = hangingStop;
integrator2.runOnGPU = runOnGPU;

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
integrator3.hangingStop = hangingStop;
integrator3.runOnGPU = runOnGPU;

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
settings.FramesToRecord = 0.5/h;
settings.InitialWindowPosition = [0, 0, 1920, 1080];
settings.VideoOut = path+string(settings.SceneName)+"/"+string(percentageImprovement);
settings.residualHistBinSize = 0.0005;

% td = simulate( {mesh2da,mesh2da2, mesh2da3}, {integrator,integrator2,integrator3}, h, settings, rigid, pcf, NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
td = simulate( {mesh2da, mesh2da3}, {integrator,integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
save(path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName)+".mat", 'td');
    % writeTDcsv(td, path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName), ["_elastic","_xpbdLayer", "_xpbdResLayer"]);
 writeTDcsv(td, path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName), ["_elastic", "_xpbdResLayer"]);  
end