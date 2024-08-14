clear;

close all;
h = 0.01; % time step

rho = 100;
nu = 0.35; % Poisson ratio: close to 0.5 for rubber
E = 2e4; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.1;
alpha1 = 0;
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

settings = SimulationSettings();
settings.SceneName = "triTest";

resetMesh = true;
scale = [1,1];
rot = 0;

angularVelocity = pi/6;
velocity = [1000,500];
% velocity = [0,0];

% mesh2d=fetchPoly2D('singleTri',resetMesh, material, scale, rot, settings);
mesh2d=fetchPoly2D('wheel.1',resetMesh, material, scale, rot, settings);
mesh2d.setRigidTransform(90,[5,5],true);
mesh2d.setRigidMotion(angularVelocity , velocity, h, true);
mesh2da = AdaptiveMesh( mesh2d);

%Doing it before or after doesn't change the error
% mesh2da.setRigidMotion(angularVelocity , velocity,h);

[R, COM, angle, inertia, rigidMasses, ri, COMDOT, angularV_out]=AdaptiveMesh.makeRigidsFromConnectivity(1, mesh2da.p, mesh2da.v, ones(mesh2da.N,1), mesh2da.mass);
disp("angular velocity diff:"+string(angularVelocity - angularV_out));
vdiff = velocity - COMDOT;
disp("velocity diff:"+string(vdiff(1))+ "," + string(vdiff(2)));