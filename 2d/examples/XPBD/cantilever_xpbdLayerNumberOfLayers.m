close all;
clear

gravity = -9.8; % gravity
h = 0.01; % time step

rho = 1;
nu = 0.37; % Poisson ratio: close to 0.5 for rubber
% E = 6e3; % Young's modulus: 0.01e9 approximate for rubber
E = 1e6;
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('cantileverP2',resetMesh, material, scale, rot, settings);

% pin the left hand side
px = mesh2d.p(1:2:end);
minx = min( px );
mesh2d.pin( find( px < minx + 0.1 ) );

offset = 3.5;
mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );
mesh2da2.plotOffset = [offset,0];
mesh2da3 = AdaptiveMesh( mesh2d );
mesh2da3.plotOffset = [offset*2,0];
mesh2da4 = AdaptiveMesh( mesh2d );
mesh2da4.plotOffset = [offset*3,0];
mesh2da5 = AdaptiveMesh( mesh2d );
mesh2da5.plotOffset = [offset*4,0];
mesh2da6 = AdaptiveMesh( mesh2d );
mesh2da6.plotOffset = [offset*5,0];
settings.CamPadding(3) = 2;
settings.CamPadding(2) = 12;
settings.campos = [0,1];

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 20; %using 200 when testing features
% iterations = 40;%using 40 for fair comparison
energyType = 0;
% computeResidual = true;
computeResidual = true;
useMassSplitting = false;
giveUpEnabled = false;
elasticCompliance = 5e-2;
runOnGPU = false;

integrator2 = xpbdResLayer2D();
integrator2.energyType = energyType;
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = 0.0;
layers = [50,0];
integrator2.setLayers(layers);
rigidIT = [1];
integrator2.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator2.computeResiduals = computeResidual;
integrator2.useMassSplitting = useMassSplitting;
integrator2.giveUpEnabled = giveUpEnabled;
integrator2.elasticCompliance = elasticCompliance;
integrator2.runOnGPU = runOnGPU;

integrator3 = xpbdResLayer2D();
% integrator3 = xpbdLayer2D();
integrator3.energyType = energyType;
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = 0.0;
layers = [66,33,0];
integrator3.setLayers(layers);
rigidIT = [1,2];
integrator3.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator3.computeResiduals = computeResidual;
integrator3.useMassSplitting = useMassSplitting;
integrator3.giveUpEnabled = giveUpEnabled;
integrator3.elasticCompliance = elasticCompliance;
integrator3.runOnGPU = runOnGPU;

integrator4 = xpbdResLayer2D();
integrator4.energyType = energyType;
integrator4.Gravity = gravity;
integrator4.boundaryCompliance = 0.0;
layers = [80,60,40,20,0];
integrator4.setLayers(layers);
rigidIT = [1,2,4,8];
integrator4.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator4.computeResiduals = computeResidual;
integrator4.useMassSplitting = useMassSplitting;
integrator4.giveUpEnabled = giveUpEnabled;
integrator4.elasticCompliance = elasticCompliance;
integrator4.runOnGPU = runOnGPU;

integrator5 = xpbdResLayer2D();
integrator5.energyType = energyType;
integrator5.Gravity = gravity;
integrator5.boundaryCompliance = 0.0;
layers = [90,80,70,60,50,40,30,20,10,0];
integrator5.setLayers(layers);
rigidIT = [2,2,2,2,2,2,2,2,2];
integrator5.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator5.computeResiduals = computeResidual;
integrator5.useMassSplitting = useMassSplitting;
integrator5.giveUpEnabled = giveUpEnabled;
integrator5.elasticCompliance = elasticCompliance;
integrator5.runOnGPU = runOnGPU;

settings.MakeVideo = 1;
settings.SceneName = 'cantilever_xpbdLayerNumberLayers';
settings.FramesToRecord  = 2/h;

settings.RigidificationEnabled = false;
settings.plotLayers = true;
if computeResidual
    settings.plotResidual = 12;%6, 9, 10
end
% settings.residualRange = [-5,-1];
settings.residualRange = [-1,0.5];
settings.camtarget = [6,0];
settings.overwriteComparisonPositionFrom1 = true;
settings.plotResidualFrames = iterations;
settings.residualComparisonLabels = ["d50","", "d33","", "d20","","d10",""];
settings.specialResidualOption = 2;
settings.PopUpPlot = true;
% settings.plotResidualRelativeComparison = 1;
% settings.DrawEDots = true;
% mesh2da.elGamma = ones(size(mesh2da.elGamma));
% mesh2da2.elGamma = ones(size(mesh2da.elGamma));

td = simulate( {mesh2da2,mesh2da3,mesh2da4, mesh2da5}, {integrator2, integrator3,integrator4,integrator5}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
save(string(settings.SceneName)+"_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
% td = simulate( {mesh2da,mesh2da2,mesh2da3}, {integrator, integrator2, integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2da,mesh2da2}, {integrator,integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2da}, {integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
save(string(settings.SceneName)+"_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td, string(settings.SceneName), settings.residualComparisonLabels);
% readTDcsv([string(settings.SceneName)+"_layer.csv",string(settings.SceneName)+"_elastic.csv"],-1,h);
% readTDcsvLog([string(settings.SceneName)+"_elastic.csv",string(settings.SceneName)+"_layer.csv"],["blue","#FF6600"],-1,h);
