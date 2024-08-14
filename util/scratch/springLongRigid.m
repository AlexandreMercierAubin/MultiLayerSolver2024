residualRigid5 = runScene(5);
residualRigid10 = runScene(10);
residualRigid15 = runScene(20);
residualRigid50 = runScene(50);
residualRigid100 = runScene(100);
residualElastic = runScene(0);
figure
plot([residualRigid5;residualRigid10;residualRigid15;residualRigid50;residualRigid100;residualElastic]')
legend(["rigid5iterations","rigid10iterations","rigid15iterations","rigid50iterations","rigid100iterations","elastic"]);
xlabel("iteration");
ylabel("|C - αTilde*lambda|");

function residualNorm = runScene(rigidEndsIterations)
    % Parameters
    mass = [1,1,1,1,1];            % Mass of the vertices
    k = [1,1,1,1];               % Spring constant
    damping = 0.5;        % Damping coefficient
    timeStep = 0.01;      % Time step for simulation
    numSteps = 1;       % Number of simulation steps
    boundaryCompliance = 0;
    
    % Initialize positions and velocities
    x = [-3,-1,1,3,5];
    v  = zeros(size(x));
    
    % Number of springs
    numSprings = numel(k);
    
    springVertices = [1,2;
                      2,3;
                      3,4;
                      4,5];
    
    restLength = [2,2,2,1.5];     % Rest length of the springs
    
    iterations = 100;
%     rigidEndsIterations = 10;
    elasticIterations = iterations - rigidEndsIterations;
    elasticSprings = [4];
    elasticVertices = unique(springVertices(elasticSprings,:));

    rigidSprings = setdiff(1:numSprings,elasticSprings);
    rigidVertices = unique(springVertices(rigidSprings,:));
    rigidBody = Rigid1D(x(rigidVertices),mass(rigidVertices),rigidVertices);
    
    vertexRigidBodyID = zeros(numel(x),1);
    vertexRigidBodyID(rigidVertices) = 1;

    rigidElasticBoundaries = intersect(elasticVertices,rigidVertices);

    % Simulation loop
    for step = 1:numSteps
        xnew = x;
        lambda = zeros(numSprings,1);
        boundaryLambda = zeros(numel(rigidElasticBoundaries),1);
        residualArray = zeros(numSprings,1);
    
        for i = 1:rigidEndsIterations
            
            [residualArray,xnew, lambda]= springConstraints1D(residualArray,elasticSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);

            for j=1:numel(rigidElasticBoundaries) %only one for now
                boundaryVert = rigidElasticBoundaries(j);
                rigidVertexPos = rigidBody.queryVertices(boundaryVert);
                rigidMass = rigidBody.mass;
                % rigidMass = mass(boundaryVert); using the elastic vertex
                % mass might lead to better convergence. This needs to be
                % tested
                [deltaLambda, deltaxElastic, rigidCOM] = rigidBoundary1D(timeStep, xnew(boundaryVert), rigidVertexPos, mass(boundaryVert), boundaryLambda(j), boundaryCompliance, rigidMass, rigidBody.com);
                rigidBody.com = rigidCOM;
                
                xnew(rigidBody.indices) = rigidBody.queryVertices(rigidBody.indices);
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

