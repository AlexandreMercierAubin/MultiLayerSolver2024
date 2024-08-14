close all
clear 

h = 0.01;
gravity = -1; % gravity

rho = 1;
nu = 0.38; % Poisson ratio: close to 0.5 for rubber
E = 1e10; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0;
alpha1 = 0.001;
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

settings = SimulationSettings();
settings.MakeVideo = 1;
% settings.CamPadding = [0 0.01 0.001 0.001];
settings.SceneName = "palm_xpbdLayer";

resetMesh = true;
scale = [50,50];
rot = 0;
mesh2d = fetchPoly2D('palm.1',resetMesh, material, scale, rot, settings);
mesh2d.setRigidTransform( 180, [ 0, 0 ] );
mesh2d.pin([1 2 3]);

mesh2da = AdaptiveMesh(mesh2d);

iterations = 40;%using 40 for fair comparison
energyType = 2;
computeResidual = true;
useMassSplitting = false;
residualType = 0; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;
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
boundaryCompliance=1e-9;

layers = [75,50,25,0];
rigidIT = [1,2,4];

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

% settings.FocusOnMeshNode = 0;
settings.campos = [0,1,0];
settings.MakeVideo = 1;
settings.plotLayers = true;
settings.PlotSkip = 1;
if computeResidual && ~ settings.plotFrameTimeHistogram
    % settings.plotResidual = 1;%6, 9, 10, 11
end
settings.residualRange = [-5,-1];
settings.overwriteComparisonPositionFrom1 = true;

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;


td = simulate( {mesh2da}, {integrator3}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());

