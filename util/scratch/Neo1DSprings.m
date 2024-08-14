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
%     rigidEndsIterations = 10;
    elasticIterations = iterations - rigidEndsIterations;
    
    % Simulation loop
    for step = 1:numSteps
        xnew = x;
        lambda = zeros(numSprings,1);
        residualArray = zeros(numSprings,1);
    
        for i = 1:rigidEndsIterations
            
            [residualArray,xnew, lambda]= springConstraints1D(residualArray,[2],springVertices,xnew, lambda, mass, restLength, k, timeStep);
    
            [deltaLambda, deltaxElastic, rigidCOM] = rigidBoundary1D(timeStep, xnew(2), xnew(1)+restLength(1), restLength(1)/2, mass(2), 0, 0, 2, xnew(1)+restLength(1)/2);
            xnew(1) = rigidCOM - restLength(1)/2;
            xnew(2) = rigidCOM + restLength(1)/2;

            [deltaLambda, deltaxElastic, rigidCOM] = rigidBoundary1D(timeStep, xnew(3), xnew(4)-restLength(3), restLength(3)/2, mass(3), 0, 0, 2, xnew(4)-restLength(3)/2);
            xnew(4) = rigidCOM + restLength(3)/2;
            xnew(3) = rigidCOM - restLength(3)/2;
        end
    
        for i = 1:elasticIterations
            [residualArray,xnew, lambda]= springConstraints1D(residualArray,1:numSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep);
        end
        v = (xnew-x)./timeStep;
        x = xnew;
    
        residualNorm = sum(residualArray(:,2:end),1);
    end
end