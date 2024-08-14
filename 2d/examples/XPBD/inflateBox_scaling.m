close all;
clear

gravity = 0; % gravity
h = 0.005; % time step

% Set damping manually
alpha0 = 0.0000; % Rayleigh factor on M
alpha1 = 0.000; % Rayleigh factor on K

%xpbd voight constraint don't match frequencies perfectly so tuning
%materials
rho = 1;
nu = 0.3; % Poisson ratio: close to 0.5 for rubber
E = 1e6; % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
material1 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

E = 3e6;
nu = 0.3;
[mu, lambda] = toLame( nu, E );
material2 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

E = 1e7;
nu = 0.4;
[mu, lambda] = toLame( nu, E );
material3 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d1 = generateCantilever2D([], material1, 0.1, scale, [3,1], 1);
mesh2d2 = generateCantilever2D([], material2, 0.065, scale, [3,1], 1);
mesh2d3 = generateCantilever2D([], material3, 0.045, scale, [3,1], 1);

%inflate
mesh2d1.p(1:2:end) = 1.1*mesh2d1.p(1:2:end);
mesh2d1.plotOffset = [-4,1];
mesh2d2.p(1:2:end) = 1.1*mesh2d2.p(1:2:end);
mesh2d2.plotOffset = [0,1];
mesh2d3.p(1:2:end) = 1.09*mesh2d3.p(1:2:end);
mesh2d3.plotOffset = [4,1];

mesh2daE = AdaptiveMesh( mesh2d1 );
mesh2da2E = AdaptiveMesh( mesh2d2 );
mesh2da3E = AdaptiveMesh( mesh2d3 );

mesh2daL = AdaptiveMesh( mesh2d1 );
mesh2da2L = AdaptiveMesh( mesh2d2 );
mesh2da3L = AdaptiveMesh( mesh2d3 );

mesh2daL.plotOffset(2)=-1;
mesh2da2L.plotOffset(2)=-1;
mesh2da3L.plotOffset(2)=-1;

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 50;
energyType = 2;
computeResidual = true;
residualType = 0; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;
if residualType == 1
    useGravityConstraints = false;
end

settings.plotFrameTimeHistogram = false;

runUntil = -1;%-1 to deactivate
giveUpEnabled = true;
hangingStop = iterations;

boundaryCompliance = 1e-9;

layers = [50,25,0];
rigidIT = [1,2];

for i =1:6
    if i <= 3
        integrators(i) = xpbdResLayer2D();
        integrators(i).iterations = [iterations];
    else
        integrators(i) = xpbdResLayer2D();
        integrators(i).setLayers(layers);
        integrators(i).iterations = [rigidIT,iterations - sum(rigidIT)];
    end
    integrators(i).energyType = energyType;
    integrators(i).Gravity = gravity;
    integrators(i).boundaryCompliance = boundaryCompliance;
    integrators(i).computeResiduals = computeResidual;
    integrators(i).giveUpEnabled = giveUpEnabled;
    integrators(i).residualType = residualType;
    integrator3.LayerOrderingType = layerType;
    integrators(i).runUntilResSmallerThan = runUntil;
    integrators(i).useGravityConstraints = useGravityConstraints;
    integrators(i).hangingStop = hangingStop;
end

settings.MakeVideo = 1;
settings.SceneName = 'inflateBox_Scaling';
settings.CamPadding = [5,5,2,2];
settings.campos = [0,0];

settings.RigidificationEnabled = false;
settings.plotLayers = true;
settings.residualRange = [-1,0.2];
% settings.DrawEDots = true;
settings.overwriteComparisonPositionFrom1 = false;
settings.FramesToRecord = 50;
settings.InitialWindowPosition = [0, 0, 1920, 1080];
settings.StopCompareToID = [1,2,3,1,2,3];
settings.StopAtImprovementPercent = 80;
if computeResidual && ~ settings.plotFrameTimeHistogram
    %settings.plotResidual = 1;%6, 9, 10, 11
end

td = simulate( {mesh2daE,mesh2da2E, mesh2da3E,mesh2daL,mesh2da2L,mesh2da3L}, {integrators(1),integrators(2),integrators(3),integrators(4),integrators(5),integrators(6)}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2daE,mesh2daL}, {integrators(1),integrators(4)}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% writeTDcsv(td, string(settings.SceneName), ["_elastic","_xpbdLayer",
% "_xpbdResLayer"]);xpb

meanComponentFindTime = [0,0,0];
meanRuntimeImprovement = [0,0,0];
for i=4:6
    meanComponentFindTime(i-3)=100*mean(td{i}.log(11,:))/mean(td{i}.log(6,:));
    meanRuntimeImprovement(i-3) = (100*mean(td{i-3}.log(6,:))/mean(td{i}.log(6,:)))-100;
end
sumTimes = [sum(td{1}.log(6,:)),sum(td{2}.log(6,:)),sum(td{3}.log(6,:)),sum(td{4}.log(6,:)),sum(td{5}.log(6,:)),sum(td{6}.log(6,:))]
meanComponentFindTime
meanRuntimeImprovement
vertices = [mesh2daE.N,mesh2da2E.N,mesh2da3E.N]
elements = [size(mesh2daE.t,1),size(mesh2da2E.t,1),size(mesh2da3E.t,1)]