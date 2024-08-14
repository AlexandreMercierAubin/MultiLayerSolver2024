close all;
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
settings.DrawDv = true;
settings.DrawApproxDv = true;
settings.MakeVideo = true;

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d=fetchPoly2D('twoTri',resetMesh, material, scale, rot, settings);

% V = [0,0;
%      0.5,0.5;
%      1,0;
%      0.5,-0.5];
% T = [2,1,3;
%     4,3,1];
% 
% mesh2d = Mesh(V,T,[],material,[1:numel(V)]',ones(numel(V),1));
% mesh2d.setRigidTransform(90,[0,0],true);

mesh2d.pin([4]);
mesh2da = AdaptiveMesh( mesh2d);
mesh2da2 = AdaptiveMesh( mesh2d);
iterations =1;
rigidIT = [1];
layers = [50];
integrator = xpbdResLayer2D();
integrator.Gravity = g; 
integrator.setLayers(layers);
integrator.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator.computeResiduals = false;
integrator.useGravityConstraints = false;

integrator2 = xpbdLayer2D();
integrator2.Gravity = g;


rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

settings.RigidificationEnabled = false;
settings.plotLayers = true;

simulate( {mesh2da,mesh2da2}, {integrator,integrator2}, h, settings, {rigid} );
% simulate( {mesh2da}, {integrator}, h, settings, {rigid} );
