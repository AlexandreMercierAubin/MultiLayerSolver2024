cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

%white
rho = 2;
nu = 0.38; 
k = 5e5;
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.0;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [1,1,1])];

%green
rho = 2;
nu = 0.39; 
k = 1e9;
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.0;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
tMaterial2 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0,1,0])];

materials = [tMaterial,tMaterial2];

scale = [1,1,1];
baseMesh = meshLoader("capsuleCoarse", [], materials, scale);

baseMesh.setRigidTransform( [0,90,0],[0,0,0]);
v = reshape(baseMesh.p,3,[]);
n1 = v( :, baseMesh.t(:,1) );
n2 = v( :, baseMesh.t(:,2) );
n3 = v( :, baseMesh.t(:,3) );
n4 = v( :, baseMesh.t(:,4) );
elCenter = 0.25 * (n1+n2+n3+n4);

dfromyaxis = elCenter(1,:);

ind = (dfromyaxis > 0);
attributes = baseMesh.materialIndex;
attributes(ind) = 2;

baseMesh.updateMaterials( attributes, materials );
baseMesh.setRigidTransform( [0,0,180],[0,-4, 2],true);
% baseMesh.setRigidTransform( [0,0,180],[0,-4, -2],true);
% baseMesh.setRigidMotion([0,pi/4,0],[0,0,0],h);

meshes = AdaptiveMesh3D(baseMesh);

gravity = -9.8;
% gravity = 9.8;

iterations = 20;

% gravity = 0;
layers = [100,50,0];
rigidIT = [1,5];
% layers = [100,0];
% rigidIT = [1];
iterationsList = [rigidIT,iterations - sum(rigidIT)];

% layers = [0];
% iterationsList = [iterations];

computeResiduals = true;
contactResistance = 1;
boundaryCompliance = 0;

energyType = 2;
computeResidual = true;
useMassSplitting = false;
residualType = 1; % 0 : constraint residual 1: eq4 residual
layerType = 0; % 0: strain rate 1: random 2:eigs
useGravityConstraints = true;

settings.plotFrameTimeHistogram = true;
if settings.plotFrameTimeHistogram
    runUntil = 1e-4;%-1 to deactivate
    giveUpEnabled = true;
    hangingStop = 20;
    % settings.StopAtImprovementPercent = 90;
else
    runUntil = 1e-4;%-1 to deactivate
    hangingStop = iterations;
    giveUpEnabled = true;
end

rigidificator = NoRigidificator3D();

integrator = xpbdLayer3D();
integrator.Gravity = gravity;
integrator.boundaryCompliance = boundaryCompliance;
integrator.contactResistance = contactResistance;
integrator.setLayers(layers);
integrator.iterations = iterationsList;
integrator.computeResiduals = computeResidual;
integrator.useMassSplitting = useMassSplitting;
integrator.giveUpEnabled = giveUpEnabled;
integrator.residualType = residualType;
integrator.LayerOrderingType = layerType;
integrator.runUntilResSmallerThan = runUntil;
integrator.useGravityConstraints = useGravityConstraints;
integrator.hangingStop = hangingStop;
integrator.airDrag = 0.005;
integrator.useRigidContacts = false;

energyModel = NeoHookean3DEnergy();

thickness = 1.3;
top = 9;
steepness = 4.25;
length = 17;

xpbdContactCompliance = 1e-4;
platform1 = PlatformContactDetector([0,0,0],[10,length,thickness],[0,0,0],0);
platform1.xpbdContactCompliance = xpbdContactCompliance;
contactFinder = {platform1};

settings.MakeVideo = 1;
settings.FramesToRecord = 2600*4;
settings.SceneName = 'pillXPBD';
% settings.PlotEDotHist = 1;
% settings.PlotSkip = 4;
settings.campos=[150,20,25]*1.5;
settings.camtarget = [0,0,3.5];
% settings.FocusOnMeshNode = [1,30];
settings.camfov = 10;
%settings.WriteOBJs = 1;
% settings.DrawRigidDOFs = 1;
settings.OBJDir = './objs/pillXPBD/';
settings.residualComparisonLabels = ["E","L"]

td = simulate3D({meshes},h,contactFinder, integrator, rigidificator, settings, energyModel);%baseMesh,
% save("pills_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
% writeTDcsv(td,"PachinkoPill",["_default", "_adaptive"]);%"_default",