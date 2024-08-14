close all;
clear

gravity = -9.8; % gravity
h = 0.01; % time step

rho = 1;
nu = 0.35; % Poisson ratio: close to 0.5 for rubber
% E = 6e3; % Young's modulus: 0.01e9 approximate for rubber
E = 1e3;
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
V = [0,0;
     0.5,0.5;
     1,0;
     0.5,-0.5];
T = [2,1,3;
    4,3,1];

mesh2d = Mesh(V,T,[],material,[1:numel(V)]',ones(numel(V),1));

% pin the left hand side
mesh2d.pin(2);

mesh2da = AdaptiveMesh( mesh2d );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 2;%using 40 for fair comparison
energyType = 2;
computeResidual = true;
useMassSplitting = false;
residualType = 0; % 0 : constraint residual 1: eq4 residual
layerType = 4; % 0: strain rate 1: random 2:eigs 3:element order 4: custom
customOrder = [2,1]';
useGravityConstraints = false;
if residualType == 1
    useGravityConstraints = false;
end

settings.plotFrameTimeHistogram = false;
hangingStop = -1;
if settings.plotFrameTimeHistogram
    runUntil = 0.1;%-1 to deactivate
    giveUpEnabled = false;
    hangingStop = 300;
else
    runUntil = -1;%-1 to deactivate
    giveUpEnabled = false;
end
% settings.residualHistBinSize = 0.0025;

% layers = [1,1];
% rigidIT = [10];

layers = [50,0];
rigidIT = [1];

integrator3 = xpbdResLayer2D();
integrator3.energyType = energyType;
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = 1e-8;
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
integrator3.customOrdering = customOrder;

settings.MakeVideo = 1;
settings.SceneName = 'cantilever_xpbdResLayerTest';
settings.CamPadding(3) = 2;
settings.campos = [0,1];

settings.RigidificationEnabled = false;
settings.plotLayers = true;
%settings.DrawGraphColors = true;

if computeResidual && ~ settings.plotFrameTimeHistogram
    % settings.plotResidual = 1;%6, 9, 10, 11
end
animation = ForceImpulseAnimationScripter();
animation.dofs = {[4*2-1;4*2]};
animation.frameNumbers = 1./h;
animation.forceImpulse = {[0,-1000]};

% settings.residualRange = [-5,-1];
settings.residualRange = [-1, 0.1];
settings.overwriteComparisonPositionFrom1 = true;
settings.skipPreProjectionResidual = false;
% settings.plotResidualRelativeComparison = 1;
% settings.DrawEDots = true;
% mesh2da.elGamma = ones(size(mesh2da.elGamma));
% mesh2da2.elGamma = ones(size(mesh2da.elGamma));
td = simulate( {mesh2da}, {integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), animation, StVenantKirchoffEnergy());
