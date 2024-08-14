close all;
clear

gravity = 0; % gravity
h = 0.01; % time step

rho = 2;
nu = 0.35; % Poisson ratio: close to 0.5 for rubber
% E = 6e3; % Young's modulus: 0.01e9 approximate for rubber
E = 1e5;
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d = AdaptiveMesh(fetchPoly2D('barP2',resetMesh, material, scale, rot, settings));
mesh2d.setRigidMotion( 1/(3*pi), [0,0],h );

mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );
mesh2da3 = AdaptiveMesh( mesh2d );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 10;%using 40 for fair comparison
energyType = 2;
computeResidual = true;
useMassSplitting = false;
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;
boundaryCompliance = 0;

settings.plotFrameTimeHistogram = false;
hangingStop = -1;
residualType= -1; % 0 : constraint residual 1: eq4 residual
if settings.plotFrameTimeHistogram
    runUntil = 0.1;%-1 to deactivate
    giveUpEnabled = false;
    hangingStop = 300;
    residualType = 1;
else
    runUntil = -1;%-1 to deactivate
    giveUpEnabled = false;
    residualType = 1;
end

if residualType == 1
    useGravityConstraints = false;
end

% settings.residualHistBinSize = 0.0025;

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
settings.SceneName = 'spin_bar_layers';
settings.CamPadding(3) = 2;
settings.campos = [0,1];
settings.FramesToRecord = 4/h;

settings.RigidificationEnabled = false;
settings.plotLayers = true;
%settings.DrawGraphColors = true;

if computeResidual && ~ settings.plotFrameTimeHistogram
    settings.plotResidual = 1;%6, 9, 10, 11
end

% settings.residualRange = [-5,-1];
settings.residualRange = [-1, 0.1];
settings.overwriteComparisonPositionFrom1 = true;
% settings.InitialWindowPosition = [0, 0, 1920, 1080];
% settings.plotResidualRelativeComparison = 1;
% settings.DrawEDots = true;
% mesh2da.elGamma = ones(size(mesh2da.elGamma));
% mesh2da2.elGamma = ones(size(mesh2da.elGamma));
settings.skipPreProjectionResidual = false;

% td = simulate( {mesh2da,mesh2da2,mesh2da3}, {integrator, integrator2, integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
td = simulate( {mesh2da,mesh2da3}, {integrator,integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2da}, {integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% save(string(settings.SceneName)+"_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td, string(settings.SceneName), ["_elastic","_xpbdLayer", "_xpbdResLayer"]);
% readTDcsv([string(settings.SceneName)+"_layer.csv",string(settings.SceneName)+"_elastic.csv"],-1,h);
% readTDcsvLog([string(settings.SceneName)+"_elastic.csv",string(settings.SceneName)+"_layer.csv"],["blue","#FF6600"],-1,h);
