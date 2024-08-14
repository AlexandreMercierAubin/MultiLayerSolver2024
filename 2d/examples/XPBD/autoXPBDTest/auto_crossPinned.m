function auto_crossPinned(path, percentageImprovement, runOnGPU)
close all;

gravity = -9.8; % gravity
h = 0.01; % time step

rho = 4;
nu = 0.35; % Poisson ratio: close to 0.5 for rubber
E = 7e5;
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();
runUntil = 0;
settings.StopAtImprovementPercent = percentageImprovement;

resetMesh = true;
scale = [1,1];
rot = 0;

% precision = 0.15;
precision = 0.4;
scale= [1,1,1];

sP = 4;
PlusBounds = [-sP,1;
    -1,1;
     -1,sP;
     1,sP;
     1,1;
     sP,1;
     sP,-1;
     1,-1;
     1,-sP;
     -1,-sP;
     -1,-1;
     -sP,-1];

fd = @(p) ddiff(ddiff(ddiff(ddiff(drectangle(p,-sP,sP,-sP,sP),drectangle(p,-sP,-1,-sP,-1)) , drectangle(p,1,sP,1,sP)) ,drectangle(p,-sP,-1,1,sP)), drectangle(p,1,sP,-sP,-1));
[p, T] = distmesh2d(fd, @huniform, precision,[-sP,-sP;sP,sP], PlusBounds);

sumpT1 = abs(p(T(:,1),1) - p(T(:,2),1)) + abs(p(T(:,1),1) - p(T(:,3),1));
isLine1 = sumpT1 <= 1e-8 ;
sumpT2 = abs(p(T(:,1),2) - p(T(:,2),2)) + abs(p(T(:,1),2) - p(T(:,3),2));
isLine2 = abs(sumpT2./3 - p(T(:,1),2)) <= 1e-8 ;

%removes elements that are in the shape of a line
T(isLine1|isLine2, :) = [];

V = [p(:,1),p(:,2)];
V(:,1) = scale(1)*V(:,1);
V(:,2) = scale(2)*V(:,2);

mesh2d = Mesh(V,T,[],material);

% pin the left hand side
px = mesh2d.p(1:2:end);
minx = min( px );
maxx = max( px );
epsilon = 0.1;
mesh2d.pin( find( px < minx + epsilon | px > maxx-epsilon) );

mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );
mesh2da3 = AdaptiveMesh( mesh2d );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

% iterations = 200; %using 200 when testing features
iterations = 40;%using 40 for fair comparison
energyType = 2;
% computeResidual = true;
computeResidual = true;
useMassSplitting = false;
residualType = 0;% 
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;
if residualType == 1
    useGravityConstraints = false;
end

settings.plotFrameTimeHistogram = true;
if settings.plotFrameTimeHistogram
    runUntil = runUntil;%-1 to deactivate
    giveUpEnabled = false;
else
    runUntil = -1;%-1 to deactivate
    giveUpEnabled = false;
end
% settings.residualHistBinSize = 0.0025;
boundaryCompliance = 1e-9;

integrator = xpbdLayer2D();
% integrator = BackwardEuler();
% integrator = xpbdLayerKinematic2D();
% integrator = xpbd2D();
integrator.Gravity = gravity;
integrator.iterations = [iterations];
% integrator.boundaryCompliance = boundaryCompliance;
% integrator.elasticCompliance = 1e-8;
integrator.computeResiduals = computeResidual;
integrator.energyType = energyType;
integrator.useMassSplitting = useMassSplitting;
integrator.residualType = residualType;
integrator.runUntilResSmallerThan = runUntil;
integrator.useGravityConstraints = useGravityConstraints;
integrator.runOnGPU = runOnGPU;

% integrator2 = xpbd2D();
integrator2 = xpbdLayer2D();
integrator2.energyType = energyType;
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = boundaryCompliance;
layers = [50,25,0];
integrator2.setLayers(layers);
rigidIT = [3,2];
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
integrator3.useGravityConstraints = useGravityConstraints;
integrator3.runOnGPU = runOnGPU;

settings.MakeVideo = 1;
settings.SceneName = 'crossPinned';
settings.CamPadding(3) = 3;
settings.CamPadding(4) = 3;
settings.CamPadding(1) = 3;
settings.campos = [0,1];

settings.RigidificationEnabled = false;
settings.plotLayers = true;
%settings.DrawGraphColors = true;

if computeResidual && ~ settings.plotFrameTimeHistogram
    settings.plotResidual = 1;%6, 9, 10, 11
end

settings.residualRange = [-1, 0.1];
settings.overwriteComparisonPositionFrom1 = true;
settings.FramesToRecord = 2./h;
settings.InitialWindowPosition = [0, 0, 1920, 1080];
settings.VideoOut = path+string(settings.SceneName)+"/"+string(percentageImprovement);

% td = simulate( {mesh2da,mesh2da2, mesh2da3}, {integrator,integrator2,integrator3}, h, settings, rigid, pcf, NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
td = simulate( {mesh2da, mesh2da3}, {integrator,integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
save(path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName)+".mat", 'td');
    % writeTDcsv(td, path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName), ["_elastic","_xpbdLayer", "_xpbdResLayer"]);
 writeTDcsv(td, path+string(settings.SceneName)+"/"+string(percentageImprovement)+string(settings.SceneName), ["_elastic", "_xpbdResLayer"]);  
end