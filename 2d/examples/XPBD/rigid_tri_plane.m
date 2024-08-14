clear;

g = -9.8;% -9.8; % gravity
h = 0.01; % time step

rho = 100;
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
mesh2d.setRigidTransform(0,[0,5],true);
mesh2da = AdaptiveMesh( mesh2d);
mesh2da2 = AdaptiveMesh( mesh2d);

integrator = xpbdLayer2D();
integrator.Gravity = g; 
integrator.iterations = 1;
integrator.useGravityConstraints = true;

integrator2 = xpbdResLayer2D();
integrator2.Gravity = g; 
layers = [100];
integrator2.setLayers(layers);
rigidIT = [1];
integrator2.iterations = [rigidIT];
integrator2.useGravityConstraints = false;


friction = 1;
pcf = PlaneContactFinder( [0.0, 1.0] , [0,4], friction );

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

settings.plotLayers = true;

simulate( {mesh2da2,mesh2da}, {integrator,integrator2}, h, settings, {rigid} ,pcf);
