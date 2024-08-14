close all;
clear

gravity = -9.8; % gravity
h = 0.008; % time step from xpbd paper

rho = 1; %mass should be 1 for all particles... so need to be enforced
nu = 0.3; %poisson ratio from xpbd paper
E = 10^5; % Young's modulus from xpbd paper
[mu, lambda] = toLame( nu, E );

% useless in this scene
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.0; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
prev = [1,2];
start = 0;
cubeSize= 1;
V = [0,cubeSize;
    0,0];
T = [];
for i = 1:10
    cubeStart = prev;
    lastVert = max(prev);
    newVert = [lastVert+1,lastVert+2];
    T = [T; 
        prev,newVert(2);
        prev(1),newVert(2),newVert(1)];
    start = start + cubeSize;
    V = [V;start,cubeSize;
        start,0];
    prev=newVert;
end
forcingMassOfVert = true(size(V(:)));
massOfVert = ones(sum(forcingMassOfVert),1);
mesh2d = Mesh( V, T, ones(size(T,1),1), material,forcingMassOfVert,massOfVert);

% pin the left hand side
mesh2d.pin( [1,2]);

mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

rigid2 = NoRigidificator();
rigid2.PreventPinnedRigidification = true;

iterations = 30;
useMassSplitting = false;
energyType = 2;
residualType = 0; % 0: constraint residual  1: xpbd paper eq 4  2: Ax-b  

integrator = xpbdLayer2D();
integrator.computeResiduals = true;
integrator.Gravity = gravity;
% integrator.separateQuicksolveGravity = false;
integrator.iterations = iterations; %iterations from xpbd paper
integrator.boundaryCompliance = 0;
%integrator.elasticCompliance = 10^(-6);%The chain particle from the paoer had 10^(-8)
integrator.layers = [0];
integrator.energyType = energyType;
integrator.useMassSplitting = useMassSplitting;
integrator.residualType = residualType; 

integrator2 = xpbdLayer2D();
integrator2.computeResiduals = true;
integrator2.Gravity = gravity;
layers = [50,25,0];
integrator2.setLayers(layers);
rigidIT = [5,5];
elasticIterations = iterations - sum(rigidIT);
assert(elasticIterations >= 0);
integrator2.iterations = [rigidIT,elasticIterations];
integrator2.boundaryCompliance = 0;
%integrator.elasticCompliance = 10^(-6);%The chain particle from the paoer had 10^(-8)
integrator2.energyType = energyType;
integrator2.useMassSplitting = useMassSplitting;
integrator2.residualType = residualType;

settings.MakeVideo = 1;
settings.SceneName = 'thin_cantilever';
% settings.DrawTimings = 1;
% settings.PrintTimings = 0;
settings.CamPadding(3) = 0;
% settings.plotTriImplicitRigidificationElastification = 1;
% settings.PlotEDotHist = 1;
% settings.RecordFramePositionInTD = true;
% settings.DrawDv = true;
% settings.DrawApproxDv = true;
% settings.StrainLimitingEnabled=true;

settings.RigidificationEnabled = false;
% settings.plotLayers = true;
% settings.DrawEDots = true;
settings.plotResidual = 1;
% settings.residualRange = [-0.5,0.5];

settings.overwriteComparisonPositionFrom1 = true;

td = simulate( {mesh2da, mesh2da2}, {integrator,integrator2}, h, settings, {rigid}, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());