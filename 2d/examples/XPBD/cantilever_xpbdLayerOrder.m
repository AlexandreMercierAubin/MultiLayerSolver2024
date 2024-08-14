close all;
clear

gravity = -9.8; % gravity
h = 0.01; % time step

rho = 1;
nu = 0.37; % Poisson ratio: close to 0.5 for rubber
% E = 5e10; %P05
E = 1e6; %P2
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

resetMesh = true;
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('cantileverP2',resetMesh, material, scale, rot, settings); %heigth 1 [-2.5,-1.5] width 3 [-1.5,1.5]

% pin the left hand side
px = mesh2d.p(1:2:end);
minx = min( px );
mesh2d.pin( find( px < minx + 0.1 ) );

offset = 3.5;
vOffset = -3;
mesh2da = AdaptiveMesh( mesh2d );
mesh2da2 = AdaptiveMesh( mesh2d );
mesh2da2.plotOffset = [offset,0];
mesh2da3 = AdaptiveMesh( mesh2d );
mesh2da3.plotOffset = [offset*2,0];
mesh2da4 = AdaptiveMesh( mesh2d );
mesh2da4.plotOffset = [0,vOffset];
mesh2da5 = AdaptiveMesh( mesh2d );
mesh2da5.plotOffset = [offset,vOffset];
mesh2da6 = AdaptiveMesh( mesh2d );
mesh2da6.plotOffset = [offset*2,vOffset];

settings.CamPadding(3) = 2;
settings.CamPadding(2) = 8;
settings.campos = [0,1];

%figuring out the stripes
p = mesh2d.getPositionFormatted()';
n1 = p( :, mesh2d.t(:,1) );
n2 = p( :, mesh2d.t(:,2) );
n3 = p( :, mesh2d.t(:,3) );
elCenter = (n1+n2+n3)/3;
vstripe1 = elCenter(2,:) < min(elCenter(2,:))+0.25;
vstripe2= (elCenter(2,:) > min(elCenter(2,:))+0.5) & (elCenter(2,:) < min(elCenter(2,:))+0.75);
elems = 1:size(mesh2da.t,1);
indsVstripes = find(vstripe1 | vstripe2);
proportionV = numel(indsVstripes)/numel(elems);
elems(indsVstripes) = [];
indsV = [indsVstripes, elems];

hstripe1 = elCenter(1,:) < min(elCenter(1,:))+0.5;
hstripe2= (elCenter(1,:) > min(elCenter(1,:))+1) & (elCenter(1,:) < min(elCenter(1,:))+1.5);
hstripe3= (elCenter(1,:) > min(elCenter(1,:))+2) & (elCenter(1,:) < min(elCenter(1,:))+2.5);
elems = 1:size(mesh2da.t,1);
indsHStripes = find(hstripe1 | hstripe2 | hstripe3);
proportionH = numel(indsHStripes)/numel(elems);
elems(indsHStripes) = [];
indsH = [indsHStripes, elems];

rigid = NoRigidificator();
rigid.PreventPinnedRigidification = true;

iterations = 20;
energyType = 2;
% computeResidual = true;
computeResidual = true;
useMassSplitting = false;
giveUpEnabled = false;
LastLayerGiveUpEnabled = false;
elasticCompliance = 5e-2;
residualType = 1;
compareResVsolver = true;
useGravityConstraints= true;
runOnGPU = false;
if residualType == 1
    useGravityConstraints = false;
end

integrator = xpbdLayer2D();
% integrator = BackwardEuler();
% integrator = xpbdLayerKinematic2D();
% integrator = xpbd2D();
integrator.Gravity = gravity;
integrator.iterations = [iterations];
% integrator.boundaryCompliance = 0;
% integrator.elasticCompliance = 1e-8;
integrator.computeResiduals = computeResidual;
integrator.energyType = energyType;
integrator.useMassSplitting = useMassSplitting;
integrator.elasticCompliance = elasticCompliance;
integrator.residualType = residualType;
integrator.useGravityConstraints = useGravityConstraints;
integrator.runOnGPU = runOnGPU;

if compareResVsolver
    integrator2 = xpbdResLayer2D();
else
    integrator2 = xpbdLayer2D();
end
integrator2.energyType = energyType;
integrator2.Gravity = gravity;
integrator2.boundaryCompliance = 0.0;
layers = [60,40,20,0];
integrator2.setLayers(layers);
rigidIT = [3,3,4];
integrator2.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator2.computeResiduals = computeResidual;
integrator2.useMassSplitting = useMassSplitting;
integrator2.giveUpEnabled = giveUpEnabled;
integrator2.elasticCompliance = elasticCompliance;
integrator2.residualType = residualType;
integrator2.LayerOrderingType=0;
integrator2.LastLayerGiveUpEnabled = LastLayerGiveUpEnabled;
integrator2.useGravityConstraints = useGravityConstraints;
integrator2.runOnGPU = runOnGPU;

if compareResVsolver
    integrator3 = xpbdResLayer2D();
else
    integrator3 = xpbdLayer2D();
end
integrator3.energyType = energyType;
integrator3.Gravity = gravity;
integrator3.boundaryCompliance = 0.0;
integrator3.setLayers(layers);
integrator3.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator3.computeResiduals = computeResidual;
integrator3.useMassSplitting = useMassSplitting;
integrator3.giveUpEnabled = giveUpEnabled;
integrator3.elasticCompliance = elasticCompliance;
integrator3.residualType = residualType;
integrator3.LayerOrderingType=2; %eigs
integrator3.LastLayerGiveUpEnabled = LastLayerGiveUpEnabled;
integrator3.useGravityConstraints = useGravityConstraints;
integrator3.runOnGPU = runOnGPU;

if compareResVsolver
    integrator4 = xpbdResLayer2D();
else
    integrator4 = xpbdLayer2D();
end
integrator4.energyType = energyType;
integrator4.Gravity = gravity;
integrator4.boundaryCompliance = 0.0;
integrator4.setLayers(layers);
integrator4.iterations = [rigidIT,iterations - sum(rigidIT)];
integrator4.computeResiduals = computeResidual;
integrator4.useMassSplitting = useMassSplitting;
integrator4.giveUpEnabled = giveUpEnabled;
integrator4.elasticCompliance = elasticCompliance;
integrator4.residualType = residualType;
integrator4.LayerOrderingType=1; %rand
integrator4.LastLayerGiveUpEnabled = LastLayerGiveUpEnabled;
integrator4.useGravityConstraints = useGravityConstraints;
integrator4.runOnGPU = runOnGPU;

if compareResVsolver
    integrator5 = xpbdResLayer2D();
else
    integrator5 = xpbdLayer2D();
end
integrator5.energyType = energyType;
integrator5.Gravity = gravity;
integrator5.boundaryCompliance = 0.0;
integrator5.setLayers([proportionV*100,0]);
integrator5.iterations = [sum(rigidIT),iterations - sum(rigidIT)];
integrator5.computeResiduals = computeResidual;
integrator5.useMassSplitting = useMassSplitting;
integrator5.giveUpEnabled = false;
integrator5.elasticCompliance = elasticCompliance;
integrator5.residualType = residualType;
integrator5.LayerOrderingType=4; %stripesvert
integrator5.LastLayerGiveUpEnabled = LastLayerGiveUpEnabled;
integrator5.customOrdering = indsV';
integrator5.useGravityConstraints = useGravityConstraints;
integrator5.runOnGPU = runOnGPU;

if compareResVsolver
    integrator6 = xpbdResLayer2D();
else
    integrator6 = xpbdLayer2D();
end
integrator6.energyType = energyType;
integrator6.Gravity = gravity;
integrator6.boundaryCompliance = 0.0;
integrator6.setLayers([proportionH*100,0]);
integrator6.iterations = [sum(rigidIT),iterations - sum(rigidIT)];
integrator6.computeResiduals = computeResidual;
integrator6.useMassSplitting = useMassSplitting;
integrator6.giveUpEnabled = false;
integrator6.elasticCompliance = elasticCompliance;
integrator6.residualType = residualType;
integrator6.LayerOrderingType=4; %stripes horizontal
integrator6.LastLayerGiveUpEnabled = LastLayerGiveUpEnabled;
integrator6.customOrdering = indsH';
integrator6.useGravityConstraints = useGravityConstraints;
integrator6.runOnGPU = runOnGPU;

settings.MakeVideo = 1;
if compareResVsolver
    settings.SceneName = 'cantilever_xpbdResLayerOrder';
else
    settings.SceneName = 'cantilever_xpbdLayerOrder';
end


settings.RigidificationEnabled = false;
settings.plotLayers = true;
if computeResidual
    settings.plotResidual = 12;%6, 9, 10, 11
end
% settings.residualRange = [-5,0.5];
% settings.residualRange = [-1,0.5];
settings.overwriteComparisonPositionFrom1 = true;
settings.plotResidualFrames = iterations;
settings.residualComparisonLabels = ["Elastic","Strain rate", "Stretch", "Random", "Vertical Stripes", "Horizontal Stripes"];
settings.FramesToRecord  = 6/h;
settings.InitialWindowPosition = [0, 0, 1920, 1080];
% settings.plotResidualRelativeComparison = 1;
% settings.DrawEDots = true;
% mesh2da.elGamma = ones(size(mesh2da.elGamma));
% mesh2da2.elGamma = ones(size(mesh2da.elGamma));
settings.skipPreProjectionResidual = true;
settings.specialResidualOption = 1;
settings.plotLayerSwitch =true;

td = simulate( {mesh2da,mesh2da2,mesh2da3,mesh2da4,mesh2da5,mesh2da6}, {integrator, integrator2, integrator3,integrator4,integrator5,integrator6}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% save("out/order/"+string(settings.SceneName)+".mat", 'td');
% td = simulate( {mesh2da,mesh2da2,mesh2da4}, {integrator, integrator2, integrator4}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2da,mesh2da2}, {integrator,integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
% td = simulate( {mesh2da}, {integrator2}, h, settings, rigid, NullContactFinder(), NullContactFinder(), NullAnimationScript(), StVenantKirchoffEnergy());
save(string(settings.SceneName)+"_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
% writeTDcsv(td, string(settings.SceneName), settings.residualComparisonLabels);
% readTDcsv([string(settings.SceneName)+"_layer.csv",string(settings.SceneName)+"_elastic.csv"],-1,h);
% readTDcsvLog([string(settings.SceneName)+"_elastic.csv",string(settings.SceneName)+"_layer.csv"],["blue","#FF6600"],-1,h);
