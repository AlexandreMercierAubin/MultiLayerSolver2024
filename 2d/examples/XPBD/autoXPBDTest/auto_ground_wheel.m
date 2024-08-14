function auto_ground_wheel(path,percentageImprovement,runOnGPU)
close all;

h = 1/100;

%TIRE
rho = 5;
nu = 0.37; % Poisson ratio: close to 0.5 for rubber
E = 5e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00001; % Rayleigh factor on M
alpha1 = 0.005; % Rayleigh factor on K

%RIM
rho2 = 5;
nu2 = 0.39; % Poisson ratio: close to 0.5 for rubber
E2 = 1e12; % Young's modulus: 0.01e9 approximate for rubber
[ mu2, lambda2 ] = toLame( nu2, E2 );

material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1,[0,0,0.5],1.2,0.90);
material2 = TriangleMaterial( rho2, mu2, lambda2, alpha0, alpha1,[1,0,0.5],1.2,0.90);
materials = [material,material2];

mesh2d = fetchPoly2D('wheel.1',true, materials);
p = mesh2d.getPositionFormatted()';
n1 = p( :, mesh2d.t(:,1) );
n2 = p( :, mesh2d.t(:,2) );
n3 = p( :, mesh2d.t(:,3) );
elCenter = (n1+n2+n3)/3;
RimDiameter = 5;
inds = find(vecnorm(elCenter',2,2)<RimDiameter);
attributes = mesh2d.materialIndex;
attributes(inds) = 2;
mesh2d.updateMaterials( attributes, materials );

mesh2d.setRigidTransform(0,[0,8],true);
mesh2d.setRigidMotion(pi/6,[0,0],h);
mesh2da = AdaptiveMesh(mesh2d);
mesh2da2 = AdaptiveMesh(mesh2d);
mesh2da3 = AdaptiveMesh(mesh2d);

settings = SimulationSettings();
runUntil = -1;
settings.StopAtImprovementPercent = percentageImprovement;

settings.DrawTimings = 0;
% settings.DrawLambdas = 1;
settings.CamPadding(3:4) = [3,1];
settings.DrawEdges = true;
settings.SceneName = "groundWheel";
settings.StrainLimitingEnabled = true;
settings.MakeVideo = 1;
settings.RigidificationEnabled = false;
settings.ElastificationEnabled = false;

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

gravity = -9.8;
layers = [100,60,30,0];
rigidIT = [1,2,4];
iterations = 30; %Note that doing too many iterations make contacts stiffer

energyType = 2;
computeResidual = true;
layerType = 0; % 0: strain rate 1: random 2:eigs

settings.plotFrameTimeHistogram = true;
giveUpEnabled = true;
hangingStop = 250;

useRigidConstacts = false;

boundaryCompliance = 1e-9;

integrator = xpbdLayer2D();
integrator.Gravity = gravity;
integrator.iterations = [iterations];
integrator.computeResiduals = computeResidual;
integrator.energyType = energyType;
integrator.runUntilResSmallerThan = runUntil;
integrator.hangingStop = hangingStop;
integrator.boundaryCompliance = boundaryCompliance;
integrator.runOnGPU = runOnGPU;

integrator2 = xpbdLayer2D();
integrator2.energyType = energyType;
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = boundaryCompliance;
integrator2.setLayers(layers);
elasticIterations = iterations - sum(rigidIT);
assert(elasticIterations >= 0);
integrator2.iterations = [rigidIT,elasticIterations];
integrator2.computeResiduals = computeResidual;
integrator2.giveUpEnabled = giveUpEnabled;
integrator2.LayerOrderingType = layerType;
integrator2.runUntilResSmallerThan = runUntil;
integrator2.hangingStop = hangingStop;
integrator2.runOnGPU = runOnGPU;
integrator2.useRigidConstacts = useRigidConstacts;

integrator3 = xpbdResLayer2D();
integrator3.energyType = energyType;
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = boundaryCompliance;
integrator3.setLayers(layers);
integrator3.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator3.computeResiduals = computeResidual;
integrator3.giveUpEnabled = giveUpEnabled;
integrator3.LayerOrderingType = layerType;
integrator3.runUntilResSmallerThan = runUntil;
integrator3.hangingStop = hangingStop;
integrator3.runOnGPU = runOnGPU;
integrator3.useRigidConstacts = useRigidConstacts;

pcf = PlaneContactFinder( [0.0, 1.0] , [0, -3], 0 );

settings.plotLayers = true;
if computeResidual && ~ settings.plotFrameTimeHistogram
    settings.plotResidual = 1;%6, 9, 10, 11
end
%settings.plotResidualFrames = iterations;
settings.overwriteComparisonPositionFrom1 = false; %overwriting means that we would include the xpbd inflation of rotational motions
settings.CamPadding = [10,5,4,3];
settings.InitialWindowPosition = [0, 0, 1920, 1080];
    settings.VideoOut = path+string(settings.SceneName)+"/"+string(percentageImprovement);
settings.FramesToRecord = 5./h;
settings.residualHistBinSize = 0.003;
settings.plotOnlyComparison = 2;

% td = simulate( {mesh2da,mesh2da2, mesh2da3}, {integrator,integrator2,integrator3}, h, settings, rigid, pcf, NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
td = simulate( {mesh2da, mesh2da3}, {integrator,integrator3}, h, settings, rigid, pcf, NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
save(path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName)+".mat", 'td');
    % writeTDcsv(td, path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName), ["_elastic","_xpbdLayer", "_xpbdResLayer"]);
 writeTDcsv(td, path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName), ["_elastic", "_xpbdResLayer"]);   
end