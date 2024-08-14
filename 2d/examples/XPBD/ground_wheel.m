clear;
close all;

h = 1/100; % time step

settings = SimulationSettings();

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

mesh2d = fetchPoly2D('wheel.1',true, materials,[1,1],[0], settings,false,[]);
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

settings.DrawTimings = 0;
% settings.DrawLambdas = 1;
settings.CamPadding(3:4) = [3,1];
settings.DrawEdges = true;
settings.SceneName = "groundWheel";
settings.MakeVideo = 1;

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

gravity = -9.8;
% gravity = 0;
layers = [100,30,0];
rigidIT = [1,2];
% layers = [100,0];
% rigidIT = [1];
iterations = 30; %Note that doing too many iterations make contacts stiffer

computeResiduals = true;

energyType = 2;
computeResidual = true;
residualType = 0; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;
runOnGPU = false;

hangingStop = -1;
settings.overwriteComparisonPositionFrom1 = false;
settings.plotFrameTimeHistogram = true;
useRigidConstacts = ~settings.overwriteComparisonPositionFrom1;

if settings.plotFrameTimeHistogram
    settings.StopAtImprovementPercent = 90;
    runUntil = 0;%-1 to deactivate
    giveUpEnabled = false;
    hangingStop = 100;
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
integrator.residualType = residualType;
integrator.runUntilResSmallerThan = runUntil;
integrator.useGravityConstraints = useGravityConstraints;
integrator.hangingStop = hangingStop;
integrator.boundaryCompliance = boundaryCompliance;
integrator.runOnGPU = runOnGPU;
integrator.useRigidConstacts = false;

integrator2 = xpbdLayer2D();
integrator2.energyType = energyType;
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = boundaryCompliance;
integrator2.setLayers(layers);
elasticIterations = iterations - sum(rigidIT);
integrator2.iterations = [rigidIT,elasticIterations];
integrator2.computeResiduals = computeResidual;
integrator2.LayerOrderingType = layerType;
integrator2.residualType = residualType;
integrator2.runUntilResSmallerThan = runUntil;
integrator2.useGravityConstraints = useGravityConstraints;
integrator2.hangingStop = hangingStop;
integrator2.runOnGPU = runOnGPU;
integrator2.useRigidConstacts = useRigidConstacts;

integrator3 = xpbdResLayer2D();
integrator3.energyType = energyType;
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = boundaryCompliance;
integrator3.setLayers(layers);
integrator3.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator3.useRigidConstacts = useRigidConstacts;

%test as purely rigid
% integrator3.setLayers([100]);
% integrator3.iterations = [iterations];

integrator3.computeResiduals = computeResidual;
integrator3.giveUpEnabled = giveUpEnabled;
integrator3.residualType = residualType;
integrator3.LayerOrderingType = layerType;
integrator3.runUntilResSmallerThan = runUntil;
integrator3.useGravityConstraints = useGravityConstraints;
integrator3.hangingStop = hangingStop;
integrator3.runOnGPU = runOnGPU;

pcf = PlaneContactFinder( [0.0, 1.0] , [0, -3], 0 );
%pcf = NullContactFinder();

settings.plotLayers = true;
if computeResidual && ~ settings.plotFrameTimeHistogram
    settings.plotResidual = 11;%6, 9, 10, 11
end
%settings.plotResidualFrames = iterations;
settings.CamPadding = [10,1,2,10];
settings.InitialWindowPosition = [0, 0, 1920, 1080];
settings.FramesToRecord = 5./h;
%settings.DrawVelocities =true;
%settings.RunRightAway = false;
% td = simulate( {mesh2da,mesh2da2, mesh2da3}, {integrator,integrator2,integrator3}, h, settings, rigid, pcf, NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
simulate( {mesh2da,mesh2da2}, {integrator,integrator3}, h, settings, rigid, pcf, NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% simulate( {mesh2da2}, {integrator3}, h, settings, rigid, pcf, NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());

% writeTDcsv(td,settings.SceneName,["elastic","xpbdLayer","xpbdResLayer"]);