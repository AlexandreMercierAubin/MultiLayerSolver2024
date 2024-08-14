[residualElastic,xelastic] = runScene([0,0]);
% residualRigid5 = runScene([5,5],xelastic);
% residualRigid10 = runScene([10,10],xelastic);
% residualRigid15 = runScene([20,20],xelastic);
% residualRigid30 = runScene([30,30],xelastic);

residualRigid5 = runScene([5,5]);
residualRigid10 = runScene([10,10]);
residualRigid15 = runScene([20,20]);
residualRigid30 = runScene([30,30]);

plot([residualRigid5;residualRigid10;residualRigid15;residualRigid30;residualElastic]')
legend(["rigid5iterations","rigid10iterations","rigid15iterations","rigid30iterations","elastic"]);
xlabel("iteration");
ylabel("|C - αTilde*lambda|");

function [residualNorm,perStepx] = runScene(rigidEndsIterations, elasticPosOverwrite)
    % Parameters
    mass = [1,1,1,1,1,1,1,1];            % Mass of the vertices
    k = [1,1,1,1,1,1,1];               % Spring constant
    damping = 0.5;        % Damping coefficient
    timeStep = 0.01;      % Time step for simulation
    numSteps = 1;       % Number of simulation steps
    boundaryCompliance = 0;
    layers = [6,2];
    
    % Initialize positions and velocities
    x = [-7,-5,-3,-1,1,3,5,7];
    
    % Number of springs
    numSprings = numel(k);
    
    springVertices = [1,2;
                      2,3;
                      3,4;
                      4,5;
                      5,6;
                      6,7;
                      7,8];
    
    restLength = [2,2,2.0,1.5,2.0,2,2];     % Rest length of the springs
    
    iterations = 100;
    elasticIterations = iterations - sum(rigidEndsIterations);
    
    prevStrain = zeros(1,numSprings);

    if nargout > 1
        perStepx = [];
    end

    % Simulation loop
    for step = 1:numSteps
        xnew = x;
        lambda = zeros(numSprings,1);
        
        residualArray = springConstraints1D(zeros(numSprings,1),1:numSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);
        residualArray = residualArray(:,2:end);

        strain = strain1D(x,springVertices,restLength);
        strainRate = (strain - prevStrain)./timeStep;

        %layers
        for layer = 1:numel(layers)
            [vertexRigidID,rigidSprings,numRigid]=clusterRigidBodies1D(strainRate,springVertices,layers(layer));
            elasticSprings = setdiff(1:numSprings,rigidSprings);
            elasticVertices = unique(springVertices(elasticSprings,:));
        
            rigidVertices = unique(springVertices(rigidSprings,:));
            rigidBodies = Rigid1D.empty;
            for i = 1:numRigid
                bodyVertices = find(vertexRigidID==i);
                rigidBodies(i) = Rigid1D(xnew(bodyVertices),mass(bodyVertices),bodyVertices);
            end
        
            rigidElasticBoundaries = intersect(elasticVertices,rigidVertices);
            boundaryLambda = zeros(numel(rigidElasticBoundaries),1);
    
            
            for i = 1:rigidEndsIterations(layer)
                
               [tmp,xnew, lambda]= springConstraints1D(residualArray,elasticSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);
        
               for j=1:numel(rigidElasticBoundaries) %only one for now
                    boundaryVert = rigidElasticBoundaries(j);
    
                    rigidBodyID = vertexRigidID(boundaryVert);
                    rigidVertexPos = rigidBodies(rigidBodyID).queryVertices(boundaryVert);
                    rigidMass = rigidBodies(rigidBodyID).mass;
                    % rigidMass = mass(boundaryVert);% using the elastic vertex
                    % mass might lead to better convergence. This needs to be
                    % tested
                    [deltaLambda, deltaxElastic, rigidCOM] = rigidBoundary1D(timeStep, xnew(boundaryVert), rigidVertexPos, mass(boundaryVert), boundaryLambda(j), boundaryCompliance, rigidMass, rigidBodies(rigidBodyID).com);
                    rigidBodies(rigidBodyID).com = rigidCOM;
                    
                    xnew(rigidBodies(rigidBodyID).indices) = rigidBodies(rigidBodyID).queryVertices(rigidBodies(rigidBodyID).indices);
                    boundaryLambda(j) = boundaryLambda(j) + deltaLambda;
               end
               residualArray = springConstraints1D(residualArray,1:numSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);
            end
        end
    
        %fully elastic layer
        for i = 1:elasticIterations
            [tmp,xnew, lambda]= springConstraints1D(residualArray,1:numSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);
            residualArray = springConstraints1D(residualArray,1:numSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);
        end

        x = xnew;

        if nargin > 1
            x = elasticPosOverwrite(step,:);
        end

        if nargout > 1
            perStepx = [perStepx;xnew];
        end
    
        residualNorm = sum(residualArray(:,2:end),1);
    end
end