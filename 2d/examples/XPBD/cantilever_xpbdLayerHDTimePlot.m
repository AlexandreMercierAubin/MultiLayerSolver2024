close all;
clear

gravity = -5; % gravity
h = 0.005; % time step

rho = 1;
nu = 0.38; % Poisson ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.0; % Rayleigh factor on K
omega0 = 6.7485; % divide by 2pi to get period
df = 0.2;  % zero is undamped, 1 is critical
% set alpha1 from alpha0
alpha0 = 0;
alpha1 = (2*df-alpha0/omega0)/omega0;

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('cantileverP05',resetMesh, material, scale, rot, settings);

% pin the left hand side
px = mesh2d.p(1:2:end);
minx = min( px );
mesh2d.pin( find( px < minx + 0.1 ) );

mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );
mesh2da3 = AdaptiveMesh( mesh2d );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 60;
computeResidual = true;
massSplitting = false;
giveupEnabled = false;

integrator = xpbdLayer2D();
integrator.Gravity = gravity;
integrator.iterations = [iterations];%6; %times 2
integrator.boundaryCompliance = 0;
% integrator.elasticCompliance = 1e-8;
integrator.computeResiduals = computeResidual;
integrator.useMassSplitting = massSplitting;

integrator2 = xpbdLayer2D();
% integrator2 = xpbdLayerKinematic2D();
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = 0;
layers = [75,50,25,0];
integrator2.setLayers(layers);
rigidIT = [5,10,15,20];
integrator2.iterations = [rigidIT,iterations - sum(rigidIT)];
% integrator2.layers = [1];
integrator2.computeResiduals = computeResidual;
integrator2.useMassSplitting = massSplitting;
integrator2.giveUpEnabled = giveupEnabled;

integrator3 = xpbdResLayer2D();
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = 0;
integrator3.setLayers(layers);
integrator3.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator3.computeResiduals = computeResidual;
integrator3.useMassSplitting = massSplitting;
integrator3.giveUpEnabled = giveupEnabled;

settings.MakeVideo = 1;
settings.SceneName = 'cantilever_layerHDTime';
% settings.DrawTimings = 1;
% settings.PrintTimings = 0;
settings.CamPadding(3) = 3;
settings.campos = [0,1.7];
% settings.plotTriImplicitRigidificationElastification = 1;
% settings.PlotEDotHist = 1;
% settings.RecordFramePositionInTD = true;
% settings.DrawDv = true;
% settings.DrawApproxDv = true;
% settings.StrainLimitingEnabled=true;

settings.RigidificationEnabled = false;
settings.plotLayers = true;
% settings.plotResidual = 1;
settings.plotResidual = 6; %6, 7, 9, 10
settings.residualRange = [-6,0];
settings.plotResidualFrames = iterations;
settings.PlotSkip = true;
% settings.DrawEDots = true;
settings.overwriteComparisonPositionFrom1 = true;


td = simulate( {mesh2da,mesh2da2}, {integrator,integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2da,mesh2da2,mesh2da3}, {integrator,integrator2,integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2da}, {integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% save(string(settings.SceneName)+"_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
% writeTDcsv(td, string(settings.SceneName), ["_elastic","_layer"]);
% readTDcsv([string(settings.SceneName)+"_layer.csv",string(settings.SceneName)+"_elastic.csv"],-1,h);
% readTDcsvLog([string(settings.SceneName)+"_elastic.csv",string(settings.SceneName)+"_layer.csv"],["blue","#FF6600"],-1,h);