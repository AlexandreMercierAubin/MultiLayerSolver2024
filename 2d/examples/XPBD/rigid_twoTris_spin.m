clear;

g = -9.8;% -9.8; % gravity
h = 0.005; % time step

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
settings.DrawDv = true;
settings.DrawApproxDv = true;
settings.MakeVideo = true;

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d=fetchPoly2D('twoTri',resetMesh, material, scale, rot, settings);
mesh2da = AdaptiveMesh( mesh2d);
mesh2da2 = AdaptiveMesh( mesh2d);
mesh2da.setRigidMotion( 0.01, [0,0] ,h);
mesh2da2.setRigidMotion( 0.01, [0,0] ,h);

integrator = xpbd2D();
erp = 0.2;
cfm = 0.1;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, erp, cfm );
integrator.Gravity = 0; 

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

rigid2 = NoRigidificator();
rigid2.PreventPinnedRigidification = true;

settings.RigidificationEnabled = false;

simulate( {mesh2da2,mesh2da}, integrator, h, settings, {rigid,rigid2} );
% simulate( {mesh2da2}, integrator, h, settings, {rigid} );
