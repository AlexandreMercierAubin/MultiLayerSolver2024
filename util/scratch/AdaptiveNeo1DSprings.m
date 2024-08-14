residualRigid5 = runScene(5);
residualRigid10 = runScene(10);
residualRigid15 = runScene(20);
residualRigid50 = runScene(50);
residualRigid100 = runScene(100);
residualElastic = runScene(0);
plot([residualRigid5;residualRigid10;residualRigid15;residualRigid50;residualRigid100;residualElastic]')
legend(["rigid5iterations","rigid10iterations","rigid15iterations","rigid50iterations","rigid100iterations","elastic"]);
xlabel("iteration");
ylabel("|C - αTilde*lambda|");

function residualNorm = runScene(rigidEndsIterations)
    % Parameters
    mass = [1.0,1,1,1];            % Mass of the vertices
    k = [1.0,1,1];               % Spring constant
    damping = 0.5;        % Damping coefficient
    timeStep = 0.01;      % Time step for simulation
    numSteps = 1;       % Number of simulation steps
    boundaryCompliance = 0;
    
    % Initialize positions and velocities
    x = [-3,-1,1,3];
    v  = zeros(size(x));
    
    % Number of springs
    numSprings = numel(k);
    
    springVertices = [1,2;
                      2,3;
                      3,4];
    
    restLength = [2.0,1.5,2.0];     % Rest length of the springs
    
    iterations = 100;
    elasticIterations = iterations - rigidEndsIterations;
    
    prevStrain = zeros(1,numSprings);

    % Simulation loop
    for step = 1:numSteps
        xnew = x;
        lambda = zeros(numSprings,1);
        
        residualArray = zeros(numSprings,1);

        strain = strain1D(x,springVertices,restLength);
        strainRate = (strain - prevStrain)./timeStep;

        [vertexRigidID,rigidSprings,numRigid]=clusterRigidBodies1D(strainRate,springVertices,2);

        elasticSprings = setdiff(1:numSprings,rigidSprings);
        elasticVertices = unique(springVertices(elasticSprings,:));
    
        rigidVertices = unique(springVertices(rigidSprings,:));
        rigidBodies = Rigid1D.empty;
        for i = 1:numRigid
            bodyVertices = find(vertexRigidID==i);
            rigidBodies(i) = Rigid1D(x(bodyVertices),mass(bodyVertices),bodyVertices);
        end
        
        vertexRigidBodyID = zeros(numel(x),1);
        vertexRigidBodyID(rigidVertices) = 1;
    
        rigidElasticBoundaries = intersect(elasticVertices,rigidVertices);
        boundaryLambda = zeros(numel(rigidElasticBoundaries),1);

        for i = 1:rigidEndsIterations
            
           [residualArray,xnew, lambda]= springConstraints1D(residualArray,elasticSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);
    
           for j=1:numel(rigidElasticBoundaries) %only one for now
                boundaryVert = rigidElasticBoundaries(j);

                rigidBodyID = vertexRigidID(boundaryVert);
                rigidVertexPos = rigidBodies(rigidBodyID).queryVertices(boundaryVert);
                rigidMass = rigidBodies(rigidBodyID).mass;
                % rigidMass = mass(boundaryVert); using the elastic vertex
                % mass might lead to better convergence. This needs to be
                % tested
                [deltaLambda, deltaxElastic, rigidCOM] = rigidBoundary1D(timeStep, xnew(boundaryVert), rigidVertexPos, mass(boundaryVert), boundaryLambda(j), boundaryCompliance, rigidMass, rigidBodies(rigidBodyID).com);
                rigidBodies(rigidBodyID).com = rigidCOM;
                
                xnew(rigidBodies(rigidBodyID).indices) = rigidBodies(rigidBodyID).queryVertices(rigidBodies(rigidBodyID).indices);
                boundaryLambda(j) = boundaryLambda(j) + deltaLambda;
            end
        end
    
        for i = 1:elasticIterations
            [residualArray,xnew, lambda]= springConstraints1D(residualArray,1:numSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);
        end
        v = (xnew-x)./timeStep;
        x = xnew;
    
        residualNorm = sum(residualArray(:,2:end),1);
    end
end