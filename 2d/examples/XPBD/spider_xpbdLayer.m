close all;
clear

gravity = -9.8; % gravity
h = 0.01; % time step

rho = 1;
nu = 0.35; % Poisson ratio: close to 0.5 for rubber
E = 1e7; % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d = AdaptiveMesh(fetchPoly2D('spider.1',resetMesh, material, scale, rot, settings));

% pin the left hand side
px = mesh2d.p(1:2:end);
py = mesh2d.p(2:2:end);
minx = min( px ) + 0.17;
maxy = max(py) - 0.5;
mesh2d.pin( find( px < minx  & py > maxy) );

mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );
mesh2da3 = AdaptiveMesh( mesh2d );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 20;
energyType = 2;
computeResidual = true;
useMassSplitting = false;
residualType = 0; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;
if residualType == 1
    useGravityConstraints = false;
end

hangingStop = -1;
% settings.residualHistBinSize = 0.0025;
boundaryCompliance = 1e-9;

giveUpEnabled = true;
runUntil = -1;

integrator3 = xpbdResLayer2D();
integrator3.energyType = energyType;
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = boundaryCompliance;
layers = [75,50,25,0];
rigidIT = [1,2,4];
integrator3.setLayers(layers);
integrator3.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator3.computeResiduals = computeResidual;
integrator3.useMassSplitting = useMassSplitting;
integrator3.giveUpEnabled = giveUpEnabled;
integrator3.residualType = residualType;
integrator3.LayerOrderingType = layerType;
integrator3.runUntilResSmallerThan = runUntil;
integrator3.useGravityConstraints = useGravityConstraints;
integrator3.hangingStop =hangingStop;

settings.MakeVideo = 1;
settings.SceneName = 'spider_xpbdLayer';
settings.CamPadding(3) = 3;
settings.campos = [0,1];

settings.RigidificationEnabled = false;
settings.plotLayers = false;
settings.DrawGraphColors = true;
% settings.plotResidual = 3;
% settings.residualRange = [-3,1];
settings.residualRange = [-1,1];
% settings.DrawEDots = true;
settings.overwriteComparisonPositionFrom1 = true;
settings.FramesToRecord = 2./h;
% settings.InitialWindowPosition = [0, 0, 1920, 1080];

% td = simulate( {mesh2da,mesh2da2, mesh2da3}, {integrator,integrator2,integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
td = simulate( {mesh2da}, {integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
writeTDcsv(td, string(settings.SceneName), ["_elastic","_xpbdLayer", "_xpbdResLayer"]);