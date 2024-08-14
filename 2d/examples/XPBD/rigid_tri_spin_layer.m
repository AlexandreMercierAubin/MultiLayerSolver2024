clear;

g = -9.8;% -9.8; % gravity
h = 0.01; % time step

rho = 10;
nu = 0.35; % Poisson ratio: close to 0.5 for rubber
E = 2e4; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.1;
alpha1 = 0;
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

settings = SimulationSettings();
% settings.RunRightAway = 0;
settings.SceneName = "triPinned";
settings.recomputeCacheAinv = true;
settings.CamPadding = [2,2,2,2];

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d=fetchPoly2D('singleTri',resetMesh, material, scale, rot, settings);
mesh2d.setRigidMotion( pi/4, [0,0],h);
mesh2da = AdaptiveMesh( mesh2d);
mesh2da2 = AdaptiveMesh( mesh2d);

computeResidual = true;
iterations = 5;

integrator = xpbdLayer2D();
integrator.Gravity = 0; 
integrator.iterations = iterations;
integrator.computeResiduals = computeResidual;

integrator2 = xpbdResLayer2D();
integrator2.Gravity = 0; 
integrator2.iterations = iterations;
layers = [100,0];
integrator2.setLayers(layers);
rigidIT = [1];
integrator2.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator2.computeResiduals = computeResidual;

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

settings.overwriteComparisonPositionFrom1 =false;
settings.plotLayers = true;
if computeResidual && ~ settings.plotFrameTimeHistogram
    settings.plotResidual = 1;%6, 9, 10, 11
end

simulate( {mesh2da2,mesh2da}, {integrator,integrator2}, h, settings, rigid );
% simulate( {mesh2da}, integrator2, h, settings, {rigid} );
