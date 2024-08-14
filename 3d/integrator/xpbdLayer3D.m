classdef xpbdLayer3D < Integrator
    % xpbd time integrator. As of now it assumes that you are using solely
    % adaptive meshes

    properties
        % none needed
        prevWarmstarts
        regularizator = 1;
        useGamma = 0;
        constraints = {};
        iterations = [10];
        %delete later, only for quick plotting of colours as an example
        layers=[0];
        sortedLayers=[0];
        layersMap = [1];
        elasticCompliance= 1e-9;
        elasticBeta = 0;
        boundaryCompliance = 1e-8;
        computeResiduals = false;
        useMassSplitting = true;
        mockBEpositions = [];
        energyType = 2; %0: STVK 1:NeoHookean 2:stvk voigt
        giveUpEnabled = false; %requires computeResiduals = true for now
        LastLayerGiveUpEnabled = false;
        giveUpThreshold = 0;
        contactResistance = 0.5;%contact damping where 0.5 is the full contact assuming two decoupled contacts are created per actual contact
        gravityCompliance = 0.001;
        residualType = 0; %0: constraint residual 1: eq4 residual 2:Ax-b
        LayerOrderingType = 0; %0: strain rate 1: random increment 2: eigs MinvK 3: element order 4:custom
        runUntilResSmallerThan = -1; %<0: stops at the given number of iterations. >0: tries to reach desired residual
        hangingStop = 1000; %if runUntil does more than this number of iterations it can be considered as hanging so we stop
        useGravityConstraints = true;
        customOrdering = [];
        alphaTildeMatrices = [];
        runOnGPU = false;
        airDrag = 0;
        useRigidContacts = true;
        StopAtImprovementPercent = 0;
        pinnedx = [];
    end
    
    methods(Static=true)

        function [xi, sumNonMonitoredTime, residualArray,runUntilResSmallerThan] = layerSolve(isElasticElement, ...
                                                                                    isRigidElement, ...
                                                                                    rigidBodySetsPerLayer, ...
                                                                                    vertexRigidBody, ...
                                                                                    isBoundaryVertex, ...
                                                                                    isRigidVertex, ...
                                                                                    xi, ...
                                                                                    lambdai, ...
                                                                                    sumNonMonitoredTime, ...
                                                                                    h , ...
                                                                                    numRigids, ...
                                                                                    residualArray, ...
                                                                                    boundaryCompliance,...
                                                                                    contactResistance,...
                                                                                    contactCompliance,...
                                                                                    iterations,...
                                                                                    layers,...
                                                                                    runUntilResSmallerThan,...
                                                                                    computeResiduals,...
                                                                                    hangingStop,...
                                                                                    giveUpEnabled,...
                                                                                    giveUpThreshold,...
                                                                                    LastLayerGiveUpEnabled,...
                                                                                    useGravityConstraints,...
                                                                                    oldp,...
                                                                                    oldOldp,...
                                                                                    contactNormals,...
                                                                                    pointColliders ,...
                                                                                    contactVertexID,...
                                                                                    mass,...
                                                                                    T,...
                                                                                    unpinnedDOFs,...
                                                                                    pinnedDOFs,...
                                                                                    DmInv,...
                                                                                    perElementXPBDalphaMatrix,...
                                                                                    perElementXPBDalphaInvMatrix,...
                                                                                    elementColor,...
                                                                                    elementDOFs,...
                                                                                    numColors,...
                                                                                    f,...
                                                                                    layersMap,...
                                                                                    gravityCompliance, ...
                                                                                    gravity, ...
                                                                                    residualv, ...
                                                                                    perElementXPBDBeta,...
                                                                                    volume,...
                                                                                    runOnGPU,...
                                                                                    useRigidContacts,...
                                                                                    StopAtImprovementPercent,...
                                                                                    pinnedx)
            errorChangeDecreased = false;
            toGPUtic = tic;
            if runOnGPU
                [residualv, isElasticElement,isRigidElement,rigidBodySetsPerLayer,vertexRigidBody,isBoundaryVertex,isRigidVertex, xi, lambdai, h , numRigids, residualArray, boundaryCompliance, contactResistance, contactCompliance, iterations, layers, runUntilResSmallerThan, computeResiduals, hangingStop, giveUpThreshold, oldp, oldOldp, contactNormals, pointColliders, contactVertexID, mass, T, unpinnedDOFs, pinnedDOFs, DmInv, perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix, elementColor, elementDOFs, numColors, f, layersMap, gravityCompliance, gravity, perElementXPBDBeta]=toGPUArray(residualv, isElasticElement, isRigidElement, rigidBodySetsPerLayer, vertexRigidBody, isBoundaryVertex, isRigidVertex, xi, lambdai, h , numRigids, residualArray, boundaryCompliance, contactResistance, contactCompliance, iterations, layers, runUntilResSmallerThan, computeResiduals, hangingStop, giveUpThreshold, oldp, oldOldp, contactNormals, pointColliders, contactVertexID, mass, T, unpinnedDOFs, pinnedDOFs, DmInv,perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix, elementColor, elementDOFs, numColors, f, layersMap, gravityCompliance, gravity, perElementXPBDBeta);
                lambdaGravity = zeros(numel(unpinnedDOFs)/3,1,'gpuArray');
                contactLambdas = zeros(numel(contactVertexID),1,'gpuArray');
            else
                lambdaGravity = zeros(numel(unpinnedDOFs)/3,1);
                contactLambdas = zeros(numel(contactVertexID),1);
            end
            sumNonMonitoredTime = sumNonMonitoredTime + toc(toGPUtic);

            LayerPropagation =-1;

            coder.gpu.nokernel();
            for layer = 1:numel(layers)
                if layer == 1
                    setPercentStop = StopAtImprovementPercent;%convoluted way to get my runUntil variable set to a percentage according to the first residual.
                else 
                    setPercentStop = 0;
                end

                xprev = xi; %not to be confused with the previous step positions oldp

                isLastLayer = layer == numel(layers);
                mappedLayerIndex = layersMap(layer);
                elasticElements = find(isElasticElement(:,mappedLayerIndex));
                rigidElements = find(isRigidElement(:,mappedLayerIndex));
                rigidIDbyElement = rigidBodySetsPerLayer(:,mappedLayerIndex);
                rigidIDbyVert = vertexRigidBody(:,mappedLayerIndex);
                isVertBoundary = isBoundaryVertex(:,mappedLayerIndex);
                rigidVerts = find(isRigidVertex(:,mappedLayerIndex));
                elasticVerts = find(~isRigidVertex(:,mappedLayerIndex));
                ElasticDOFs = reshape([elasticVerts*3-2;elasticVerts*3-1;elasticVerts*3],[],1);

                isComponentPinned = false(numRigids(mappedLayerIndex), 1);
                [R, COM, inertia0, rigidMasses, ri, ~, COMDOT,angularv]=AdaptiveMesh3D.makeRigidsFromConnectivity(numRigids(mappedLayerIndex), xi, residualv, rigidVerts, isComponentPinned, rigidIDbyVert, rigidIDbyElement, ElasticDOFs, mass, unpinnedDOFs);
                inertia = inertia0;
                if numel(rigidMasses) == 0 && numel(rigidElements)>0
                    disp("failure to build any rigid body");
                end

                %step rigids and released particles
                %release all elastic residual
                xi(ElasticDOFs) = xi(ElasticDOFs) + h*residualv(ElasticDOFs);
                residualv(ElasticDOFs) = 0;

                steppedVelocityR = R;
                if (layer <= LayerPropagation) || (LayerPropagation == -1 && numel(rigidElements)>0)
                    % step rigids
                    [steppedVelocityR, COM, inertia] = RigidBody3D.rigidSymplecticImperativeStep(COM, R, h, COMDOT, angularv, inertia0);

                    %update pos wrt stepped rigids
                    xi = xpbdLayer3D.rotateAllVerts(rigidIDbyVert, ri, steppedVelocityR, COM, xi);
                    
                    isRigidVert = rigidIDbyVert > 0;
                    rigidInds = find(isRigidVert);
                    rigidDOFs = reshape([rigidInds*3-2, rigidInds*3-1, rigidInds*3]',[],1);

                    rigidResidualv = zeros(size(xi,1),1);
                    rigidResidualv(rigidDOFs) = BDF.BDF1(xi(rigidDOFs), xprev(rigidDOFs),h);

                    rotatedResv = xpbdLayer3D.rotateAllVerts(rigidIDbyVert,formatPositions3D(residualv),steppedVelocityR,zeros(1,3,size(steppedVelocityR,3)),zeros(size(xi)));

                    residualv(rigidDOFs) = rotatedResv(rigidDOFs) - rigidResidualv(rigidDOFs);
                    
                end
                normresv = norm(residualv);

                [xi, lambdai,  lambdaGravity, residualArray,sumNonMonitoredTime,errorChangeDecreased, stopSolve, projectedR, ~, ~, ~,runUntilResSmallerThan] = xpbdLayer3D.projectGS(xi,...
                                                                                                                                        lambdai, ...
                                                                                                                                        lambdaGravity, ...
                                                                                                                                        h, ...
                                                                                                                                        elasticElements, ...
                                                                                                                                        isVertBoundary, ...
                                                                                                                                        residualArray, ...
                                                                                                                                        sumNonMonitoredTime, ...
                                                                                                                                        errorChangeDecreased, ...
                                                                                                                                        contactLambdas, ...
                                                                                                                                        boundaryCompliance,...
                                                                                                                                        contactResistance,...
                                                                                                                                        contactCompliance,...
                                                                                                                                        iterations(layer),...
                                                                                                                                        runUntilResSmallerThan,...
                                                                                                                                        computeResiduals,...
                                                                                                                                        hangingStop,...
                                                                                                                                        giveUpEnabled,...
                                                                                                                                        giveUpThreshold,...
                                                                                                                                        isLastLayer,...
                                                                                                                                        LastLayerGiveUpEnabled,...
                                                                                                                                        useGravityConstraints,...
                                                                                                                                        oldp,...
                                                                                                                                        oldOldp,...
                                                                                                                                        contactNormals,...
                                                                                                                                        pointColliders ,...
                                                                                                                                        contactVertexID,...
                                                                                                                                        mass,...
                                                                                                                                        T,...
                                                                                                                                        unpinnedDOFs,...
                                                                                                                                        pinnedDOFs,...
                                                                                                                                        rigidIDbyVert,...
                                                                                                                                        DmInv,...
                                                                                                                                        perElementXPBDalphaMatrix,...
                                                                                                                                        perElementXPBDalphaInvMatrix,...
                                                                                                                                        elementColor,...
                                                                                                                                        elementDOFs,...
                                                                                                                                        numColors,...
                                                                                                                                        f, ...
                                                                                                                                        steppedVelocityR, ...
                                                                                                                                        COM, ...
                                                                                                                                        inertia, ...
                                                                                                                                        inertia0, ...
                                                                                                                                        rigidMasses, ...
                                                                                                                                        ri,...
                                                                                                                                        gravityCompliance, ...
                                                                                                                                        gravity, ...
                                                                                                                                        perElementXPBDBeta,...
                                                                                                                                        volume,...
                                                                                                                                        runOnGPU,...
                                                                                                                                        normresv,...
                                                                                                                                        useRigidContacts,...
                                                                                                                                        setPercentStop,...
                                                                                                                                        pinnedx);
                if (layer <= LayerPropagation) || (LayerPropagation == -1 && numel(rigidElements)>0)
                    constraintedR = cat(3,projectedR);
                    deltaR = pagemtimes(constraintedR,'none',steppedVelocityR,'transpose');
                    residualv = xpbdLayer3D.rotateResidual(rigidIDbyVert, residualv, deltaR);
                end

                if stopSolve
                    break;
                end
            end
            if runOnGPU
                fromGPU = tic;
                [xi, sumNonMonitoredTime, residualArray]=toCPUarray(xi, sumNonMonitoredTime, residualArray);
                sumNonMonitoredTime = sumNonMonitoredTime + toc(fromGPU);
            end
        end

        function residualv = rotateResidual(rigidIDbyVert, residualv, Rstack) 
            %Rstack is a 3D array of rotation matrices 3by3byNrigids
            rigidVerts = find(rigidIDbyVert ~= 0);
            resvFormatted = pagePositions3D(residualv);
            resvCropped = resvFormatted(:,:,rigidVerts);
            newResv = pagemtimes(Rstack(:,:,rigidIDbyVert(rigidVerts)),'none',resvCropped,'transpose');
            residualv([rigidVerts*3-2,rigidVerts*3-1,rigidVerts*3]) = permute(newResv,[3,1,2]);
        end

        function [xi, lambdai,  lambdaGravity, residualArray,sumNonMonitoredTime,errorChangeDecreased, stopSolve, R, COM, inertia, rigidMasses,runUntilResSmallerThan] = projectGS(xi, ...
                                                                                                                                                        lambdai, ...
                                                                                                                                                        lambdaGravity, ...
                                                                                                                                                        h, ...
                                                                                                                                                        elasticElements, ...
                                                                                                                                                        isVertBoundary, ...
                                                                                                                                                        residualArray, ...
                                                                                                                                                        sumNonMonitoredTime, ...
                                                                                                                                                        errorChangeDecreased, ...
                                                                                                                                                        contactLambdas, ...
                                                                                                                                                        boundaryCompliance,...
                                                                                                                                                        contactResistance,...
                                                                                                                                                        contactCompliance,...
                                                                                                                                                        iterations,...
                                                                                                                                                        runUntilResSmallerThan,...
                                                                                                                                                        computeResiduals,...
                                                                                                                                                        hangingStop,...
                                                                                                                                                        giveUpEnabled,...
                                                                                                                                                        giveUpThreshold,...
                                                                                                                                                        isLastLayer,...
                                                                                                                                                        LastLayerGiveUpEnabled,...
                                                                                                                                                        useGravityConstraints, ...
                                                                                                                                                        oldp,...
                                                                                                                                                        oldOldp,...
                                                                                                                                                        contactNormals,...
                                                                                                                                                        pointColliders ,...
                                                                                                                                                        contactVertexID,...
                                                                                                                                                        mass,...
                                                                                                                                                        T,...
                                                                                                                                                        unpinnedDOFs,...
                                                                                                                                                        pinnedDOFs,...
                                                                                                                                                        rigidIDbyVert,...
                                                                                                                                                        DmInv,...
                                                                                                                                                        perElementXPBDalphaMatrix,...
                                                                                                                                                        perElementXPBDalphaInvMatrix,...
                                                                                                                                                        elementColor,...
                                                                                                                                                        elementDOFs,...
                                                                                                                                                        numColors,...
                                                                                                                                                        f, ...
                                                                                                                                                        R, ...
                                                                                                                                                        COM, ...
                                                                                                                                                        inertia, ...
                                                                                                                                                        inertia0,...
                                                                                                                                                        rigidMasses, ...
                                                                                                                                                        ri,...
                                                                                                                                                        gravityCompliance, ...
                                                                                                                                                        gravity, ...
                                                                                                                                                        perElementXPBDBeta,...
                                                                                                                                                        volume,...
                                                                                                                                                        runOnGPU,...
                                                                                                                                                        normresv,...
                                                                                                                                                        useRigidContacts,...
                                                                                                                                                        setPercentStop,...
                                                                                                                                                        pinnedx)
    
            
            boundaryVertToSolve = find(isVertBoundary);

            toGPUtimer = tic;
            if runOnGPU
                lambdaBoundary = zeros(numel(boundaryVertToSolve),1,"gpuArray");
                boundaryVertToSolve = gpuArray(boundaryVertToSolve);
            else
                lambdaBoundary = zeros(numel(boundaryVertToSolve),1);
            end
            sumNonMonitoredTime = sumNonMonitoredTime + toc(toGPUtimer);
            
            %GS solve ----
            stopSolve = false;
            it = 1;

            coder.gpu.nokernel();
            while it <= iterations || (~stopSolve && runUntilResSmallerThan >= 0 && isLastLayer)
                it = it + 1;

                [lambdai,xi] = xpbdLayer3D.stvkVoigtConstraints(numColors, elementDOFs, elementColor, DmInv, perElementXPBDalphaMatrix./(h*h), mass, volume, lambdai,xi,elasticElements, perElementXPBDBeta.*(h*h), oldp, h);


                %rigid-elastic glue distance constraints 
                [lambdaBoundary, xi, R, COM, inertia, rigidMasses] = xpbdLayer3D.rigidElasticBoundaryAccumulationConstraints(boundaryCompliance, mass, rigidIDbyVert, R, COM, inertia, inertia0, rigidMasses, ri, xi, boundaryVertToSolve, lambdaBoundary, h); 

                %handling gravity
                if useGravityConstraints
                    % [lambdaGravity, xi,~] = xpbdLayer3D.gravityConstraints(mass, unpinnedDOFs, lambdaGravity, xi, h, oldp, gravityCompliance, gravity);
                end

                %to handle contacts as hard constraints
                [contactLambdas, xi, R, COM, inertia, rigidMasses] = xpbdLayer3D.contactConstraints(contactNormals, pointColliders, contactVertexID, mass, rigidIDbyVert, contactLambdas, xi, h, contactResistance, contactCompliance, R, COM, inertia, inertia0, rigidMasses, ri, useRigidContacts);

                %handling pinned constraints
                xi(pinnedDOFs) = pinnedx;

                ticNonMonitoredTime = tic;
                [endLayer, residualArray, errorChangeDecreased, ~, stopSolve,setPercentStop,runUntilResSmallerThan] = xpbdLayer3D.simpleMonitorError(h, lambdai, xi, oldp, oldOldp, mass, T, f, pinnedDOFs, elementDOFs, DmInv, perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix, residualArray, errorChangeDecreased, computeResiduals, runUntilResSmallerThan, hangingStop, giveUpEnabled, giveUpThreshold, isLastLayer, LastLayerGiveUpEnabled, lambdaGravity, gravityCompliance, gravity, unpinnedDOFs, contactNormals, pointColliders, contactVertexID, contactLambdas, contactResistance, contactCompliance, ri, perElementXPBDBeta, volume, normresv,setPercentStop);
                sumNonMonitoredTime =  toc(ticNonMonitoredTime) + sumNonMonitoredTime;
                if endLayer || stopSolve
                    break;
                end
            end

            %because no rigid body share dofs, we can vectorize
            xi = xpbdLayer3D.rotateAllVerts(rigidIDbyVert, ri, R, COM, xi);
        end

        function [lambdai,xi] = stvkVoigtConstraints(numColors, elementDOFs, elementColor, DmInv, alphaTilde, mass, volume, lambdai,xi,elasticElements, betaTilde, oldp, h)
            if numel(elasticElements) == 0
                return;
            end
            graphColors = elementColor(elasticElements);

            minv = 1./mass;

             %https://www.mathworks.com/help/gpucoder/ref/coder.gpu.nokernel.html
            coder.gpu.nokernel();
            for c = 1 : numColors
                elementsOfColor = elasticElements(graphColors==c);                
                colorDOFs = elementDOFs(elementsOfColor,:);

                xiOfColor = reshape(xi(colorDOFs)',12,1,[]);
                xiOldOfColor = reshape(oldp(colorDOFs)',12,1,[]);
                minvOfColor = reshape(minv(colorDOFs)',12,1,[]);
                lambdaOfColor = reshape(lambdai(elementsOfColor,:)',6,1,[]);
                alphaTildeOfColor = alphaTilde(:,:,elementsOfColor);
                betaTildeOfColor = betaTilde(:,:,elementsOfColor);
                DmInvOfColor = DmInv(:,:,elementsOfColor);
                volumeOfColor = reshape(volume(elementsOfColor),1,1,[]);

                [deltaxi, ~, deltaLambdai, ~] = matStvkVoigtConstraint3D( xiOfColor, xiOldOfColor, DmInvOfColor, alphaTildeOfColor, lambdaOfColor, minvOfColor, volumeOfColor, betaTildeOfColor, h);
                lambdaOfColor = lambdaOfColor + deltaLambdai;
                xiOfColor = xiOfColor + deltaxi;

                %This overwritting is only ok because the colors are
                %hitting independant dofs
                lambdai(elementsOfColor,:)=reshape(lambdaOfColor,6,[])';
                xi(colorDOFs) = reshape(xiOfColor,12,[])';
            end
        end

        function [lambdaBoundary, xi, R, COM, inertia, rigidMasses] = rigidElasticBoundaryAccumulationConstraints(boundaryCompliance, mass, rigidIDbyVert, R, COM, inertia, inertia0, rigidMasses, ri, xi, boundaryVertToSolve, lambdaBoundary, h)%#codegen
            assert(size(COM,2)==3);
            bodyIDs = rigidIDbyVert(boundaryVertToSolve);
            boundaryElasticDOFs = reshape([boundaryVertToSolve*3-2,boundaryVertToSolve*3-1,boundaryVertToSolve*3]',1,3,[]);
            elasticxi = xi(boundaryElasticDOFs);
            massInvElastic = reshape(1./mass(boundaryVertToSolve*3),1,1,[]);

            riCropped = permute(ri(boundaryVertToSolve,:),[3,2,1]);
            rigidxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,bodyIDs), COM(:,:,bodyIDs));
            rigidInertia = inertia(:,:,bodyIDs);

            r = rigidxi-COM(:,:,bodyIDs);
            boundaryVertsCount = accumarray(bodyIDs,ones(size(bodyIDs)));
            splitUpdateByVerts = reshape(boundaryVertsCount(bodyIDs),1,1,[]);

            [deltaLambda, ~, constraintImpulse, deltaOmega] = rigidBoundaryConstraint.project(h, elasticxi, rigidxi, r, massInvElastic, reshape(lambdaBoundary,1,1,[]), boundaryCompliance, reshape(rigidMasses(bodyIDs),1,1,[]), rigidInertia, splitUpdateByVerts);

            lambdaBoundary = lambdaBoundary + reshape(deltaLambda,[],1,1);
            numRigids = size(R,3);
            if(numel(lambdaBoundary)>0)
                deltaCOMxAccumulator = reshape(accumarray(bodyIDs,reshape(constraintImpulse(:,1,:),[],1,1),[numRigids,1]),1,1,[]);
                deltaCOMyAccumulator = reshape(accumarray(bodyIDs,reshape(constraintImpulse(:,2,:),[],1,1),[numRigids,1]),1,1,[]);
                deltaCOMzAccumulator = reshape(accumarray(bodyIDs,reshape(constraintImpulse(:,3,:),[],1,1),[numRigids,1]),1,1,[]);
                deltaCOMperRigid = cat(2,deltaCOMxAccumulator,deltaCOMyAccumulator,deltaCOMzAccumulator);
                COM = COM + deltaCOMperRigid;

                omegax = reshape(accumarray(bodyIDs,reshape(deltaOmega(:,1,:),[],1,1),[numRigids,1]),1,1,[]);
                omegay = reshape(accumarray(bodyIDs,reshape(deltaOmega(:,2,:),[],1,1),[numRigids,1]),1,1,[]);
                omegaz = reshape(accumarray(bodyIDs,reshape(deltaOmega(:,3,:),[],1,1),[numRigids,1]),1,1,[]);
                constraintVelocities = pagetranspose(pagemldivide(inertia,cat(1,omegax,omegay,omegaz)));
                deltaR = expRodrigues(constraintVelocities,h);
                R = pagemtimes(deltaR,R);

                inertia = pagemtimes(R,pagemtimes(inertia0,'none',R,'transpose'));

                riCropped = permute(ri(boundaryVertToSolve,:),[3,2,1]);
                rigidxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,bodyIDs),COM(:,:,bodyIDs));
                xi(boundaryElasticDOFs) = rigidxi; %constraint is exactly enforced to preserve the rigid boundary
            end
        end

        function [contactLambdas, xi, R, COM, inertia, rigidMasses] = contactConstraints(contactNormals, pointColliders, vertexIDs, mass, rigidIDbyVert, contactLambdas, xi, h, contactResistance, contactCompliance, R, COM, inertia, inertia0, rigidMasses, ri, useRigidContacts)%#codegen
            %I probably want to paralellize these... should split between
            %rigid and elastic, then solve the rigid ones sequentially when
            %they belong to the same rigid body
            minv = 1./mass;

            isElasticContact = rigidIDbyVert(vertexIDs) == 0;
            elasticContactVerts = vertexIDs(isElasticContact);
            elasticContactPoints = pagePositions3D(pointColliders(isElasticContact,:));
            elasticContactNormals = pagePositions3D(contactNormals(isElasticContact, :));
            elasticContactDOFs = reshape([elasticContactVerts*3-2, elasticContactVerts*3-1, elasticContactVerts*3]',1,3,[]); 
            formattedxi = reshape(xi(elasticContactDOFs),1,3,[]);
            elasticminv = reshape(minv(elasticContactDOFs) ,1,3,[]);

            [dx, deltalambda,~]=penaltyConstraint.projectElastic(formattedxi, reshape(contactLambdas(isElasticContact),1,1,[]), elasticminv, elasticContactPoints, elasticContactNormals, contactResistance, reshape(contactCompliance(isElasticContact),1,1,[]), h);
            xi(elasticContactDOFs) = formattedxi + dx;
            contactLambdas(isElasticContact) = contactLambdas(isElasticContact) + deltalambda(:);       

            isRigidContact = ~isElasticContact;
            rigidContactVerts = vertexIDs(isRigidContact);
            if numel(rigidContactVerts) <1 || ~useRigidContacts
                return;
            end
            bodyIDs = rigidIDbyVert(rigidContactVerts);
            rigidParticleDOFs = reshape([rigidContactVerts*3-2, rigidContactVerts*3-1,rigidContactVerts*3]',1,3,[]);
            rigidParticlexi = reshape(xi(rigidParticleDOFs),1,3,[]);

            riCropped = permute(ri(rigidContactVerts,:),[3,2,1]);
            rigidxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,bodyIDs), COM(:,:,bodyIDs));
            rigidInertia = inertia(:,:,bodyIDs);

            r = rigidxi-COM(:,:,bodyIDs);
            boundaryVertsCount = accumarray(bodyIDs,ones(size(bodyIDs)),[size(R,3),1]);
            splitUpdateByVerts = reshape(boundaryVertsCount(bodyIDs),1,1,[]);
            rigidContactPoint = reshape(pointColliders(isRigidContact,:)',1,3,[]);
            rigidContactNormal = reshape(contactNormals(isRigidContact,:)',1,3,[]);
            rigidParticleminv = reshape(minv(rigidParticleDOFs) ,1,3,[]);

            rigidContactsScaleFactor = 2;
            [deltaLambda, deltaCOM, deltaOmega]=penaltyConstraint.projectRigid(h, rigidParticlexi, rigidContactPoint, rigidContactNormal, contactResistance, r, reshape(contactLambdas(isRigidContact),1,1,[]), reshape(rigidContactsScaleFactor*contactCompliance(isRigidContact),1,1,[]),  reshape(rigidMasses(bodyIDs),1,1,[]), rigidParticleminv(:,1,:), rigidInertia,splitUpdateByVerts);

            numRigids = size(COM,3);
            contactLambdas(isRigidContact) = contactLambdas(isRigidContact) + reshape(deltaLambda,[],1,1);
            if(numel(contactLambdas(isRigidContact))>0)
                deltaCOMxAccumulator = reshape(accumarray(bodyIDs,reshape(deltaCOM(:,1,:),[],1,1),[numRigids,1]),1,1,[]);
                deltaCOMyAccumulator = reshape(accumarray(bodyIDs,reshape(deltaCOM(:,2,:),[],1,1),[numRigids,1]),1,1,[]);
                deltaCOMzAccumulator = reshape(accumarray(bodyIDs,reshape(deltaCOM(:,3,:),[],1,1),[numRigids,1]),1,1,[]);
                deltaCOMAccumulator = cat(2,deltaCOMxAccumulator,deltaCOMyAccumulator,deltaCOMzAccumulator);
                COM = COM + deltaCOMAccumulator;

                omegax = reshape(accumarray(bodyIDs,reshape(deltaOmega(:,1,:),[],1,1),[numRigids,1]),1,1,[]);
                omegay = reshape(accumarray(bodyIDs,reshape(deltaOmega(:,2,:),[],1,1),[numRigids,1]),1,1,[]);
                omegaz = reshape(accumarray(bodyIDs,reshape(deltaOmega(:,3,:),[],1,1),[numRigids,1]),1,1,[]);
                constraintVelocities = pagetranspose(pagemldivide(inertia,cat(1,omegax,omegay,omegaz)));
                deltaR = expRodrigues(constraintVelocities,h);
                R = pagemtimes(deltaR,R);

                inertia = pagemtimes(R,pagemtimes(inertia0,'none',R,'transpose'));

                riCropped = permute(ri(rigidContactVerts,:),[3,2,1]);
                rigidxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,bodyIDs), COM(:,:,bodyIDs));
                xi(rigidParticleDOFs) = rigidxi; %constraint is exactly enforced to preserve the rigid boundary
            end
        end

        function force = scriptedForceInjection(mesh3D, animationScripts, frame, h)
            force = mesh3D.f;
            for animID = 1:numel(animationScripts)
                force = animationScripts{animID}.scriptForceAnimation(force,frame,h);
            end
        end

        function [lambdaGravity, xi, constraintC] = gravityConstraints(mass, unpinnedDOFs, lambdaGravity, xi, h, prevx, gravityCompliance, gravity) %#codegen

            dofs = formatPositions3D(unpinnedDOFs);
            verticalDOFs = dofs(:,3);
            parallelxi = xi(verticalDOFs);
            parallelxold = prevx(verticalDOFs);
            parallelMass = mass(verticalDOFs);

            coder.gpu.kernelfun();
            [deltax, deltaLambda, constraintC] = gravityConstraint.project(parallelxi, parallelxold, lambdaGravity, parallelMass, gravityCompliance , gravity ,h);
            lambdaGravity = lambdaGravity + deltaLambda;
            xi(verticalDOFs) = parallelxi + deltax;           
        end

        function xi = rotateAllVerts(rigidIDbyVert, ri, R, COM, xi)
            rigidVerts = find(rigidIDbyVert ~= 0);
            riCropped = permute(ri(rigidVerts,:),[3,2,1]);
            queriedxi = pagetranspose(RigidBody3D.pageRelativeVertexPosition(riCropped,R(:,:,rigidIDbyVert(rigidVerts)), COM(:,:,rigidIDbyVert(rigidVerts))));
            xi([rigidVerts*3-2,rigidVerts*3-1,rigidVerts*3]) = permute(queriedxi,[3,1,2]);
        end

        function [endLayer, residualArray, errorChangeDecreased, newResidual, stopSolve,setPercentStop,runUntilResSmallerThan] = simpleMonitorError(h, lambdai, xi, oldp, oldOldp, mass, T, f, pinnedDOFs, elementDOFs, DmInv,perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix, residualArray, errorChangeDecreased, computeResiduals, runUntilResSmallerThan, hangingStop, giveUpEnabled, giveUpThreshold, isLastLayer, LastLayerGiveUpEnabled, lambdaGravity, gravityCompliance, gravity, unpinnedDOFs, contactNormals, pointColliders, contactVertexIDs, contactLambdas, contactResistance, contactCompliance, ri, perElementXPBDBeta, volume, normresv,setPercentStop)
            endLayer = false; %move on to next layer
            stopSolve = false; %complete stop of solve
            newResidual = 0;
            if computeResiduals
                newRes = xpbdLayer3D.simpleEvaluateAllConstraintResiduals(xi, oldp, oldOldp, mass, T, f, pinnedDOFs, h, lambdai, lambdaGravity, gravityCompliance, gravity, unpinnedDOFs, elementDOFs, DmInv, perElementXPBDalphaMatrix, contactNormals, pointColliders, contactVertexIDs, contactLambdas, contactResistance, contactCompliance, ri, perElementXPBDBeta, volume);
                residualArray = [residualArray,newRes]; 

                lastResidual = norm(residualArray(:,end-1));
                newResidual = norm(newRes);
                errorChange = newResidual-lastResidual;

                if setPercentStop > 0
                    runUntilResSmallerThan = max((1-setPercentStop/100.0)*newResidual,1e-8);
                    setPercentStop = 0;
                end
                if errorChange < 0 %made this mechanism because sometimes the residual increases at first even in an fully elastic setup... because of this it would make the layer solver give up on layers prematurely
                    errorChangeDecreased = true;
                end

                if giveUpEnabled   && errorChange > giveUpThreshold && errorChangeDecreased && (~isLastLayer || LastLayerGiveUpEnabled) %&& numel(rigidElements) > 0 && it > 1
                    endLayer = true;
                end

                stopHanging = size(residualArray,2) > hangingStop &&  hangingStop >= 0;
                if stopHanging || (runUntilResSmallerThan >= 0 && (runUntilResSmallerThan >= newResidual && (isLastLayer||normresv < runUntilResSmallerThan)) )
                    stopSolve = true;
                end
            end
        end

        function [residualVector] = simpleEvaluateAllConstraintResiduals(xi, oldp, oldOldp, mass, T, f, pinnedDOFs, h, lambdai, lambdaGravity, gravityCompliance, gravity, unpinnedDOFs, elementDOFs, DmInv, perElementXPBDalphaMatrix, contactNormals, pointColliders, contactVertexIDs, contactLambdas, contactResistance, contactCompliance, ri, perElementXPBDBeta, volume)
            h2 = h*h;
            alphaTilde = perElementXPBDalphaMatrix./h2;
            betaTilde = perElementXPBDBeta .*h2;
            minv = 1./mass;

            coder.gpu.kernelfun();              
            %here everything is considered the same color so we can
            %evaluate in parallel
            xiOfColor = reshape(xi(elementDOFs)',12,1,[]);
            xiOldOfColor = reshape(oldp(elementDOFs)',12,1,[]);
            minvOfColor = reshape(minv(elementDOFs)',12,1,[]);
            elementsOfColor = 1:size(T,1);
            lambdaOfColor = reshape(lambdai(elementsOfColor,:)',6,1,[]);
            alphaTildeOfColor = alphaTilde(:,:,elementsOfColor);
            betaTildeOfColor = betaTilde(:,:,elementsOfColor);
            DmInvOfColor = DmInv(:,:,elementsOfColor);
            volumeOfColor = reshape(volume(elementsOfColor),1,1,[]);

            [~, ~, ~, stvkC] = matStvkVoigtConstraint3D( xiOfColor, xiOldOfColor, DmInvOfColor, alphaTildeOfColor, lambdaOfColor, minvOfColor, volumeOfColor, betaTildeOfColor, h);
            correctionTerm = abs(pagemtimes(alphaTilde,lambdaOfColor));
            absStvk = abs(stvkC);
            constraintResiduals = absStvk - min(absStvk,correctionTerm);


            %here analyses the constraint residual as purely elastic
            %contacts
            elasticContactPoints = pagePositions3D(pointColliders);
            elasticContactNormals = pagePositions3D(contactNormals);
            elasticContactDOFs = pagePositions3D([contactVertexIDs*3-2, contactVertexIDs*3-1, contactVertexIDs*3]); 
            formattedxi = reshape(xi(elasticContactDOFs),1,3,[]);
            elasticminv = reshape(minv(elasticContactDOFs) ,1,3,[]);

            [~,~,contactC]=penaltyConstraint.projectElastic(formattedxi, reshape(contactLambdas,1,1,[]), elasticminv, elasticContactPoints, elasticContactNormals, contactResistance, reshape(contactCompliance,1,1,[]), h);

            absContact = abs(contactC(:));
            contactCorrection = abs((contactCompliance/h2).*contactLambdas(:));
            constraintResidualsContacts = absContact - min(absContact,contactCorrection);

            residualVector = norm([constraintResiduals(:);constraintResidualsContacts(:)],2);
        end
    end
    
    methods
        function obj = xpbdLayer3D()
            obj@Integrator();
            obj.Name = 'Backward Euler for 3D';
        end

        function setLayers(obj, layersArray)
            obj.layers = layersArray;
            [obj.sortedLayers, ia, ic] = unique(layersArray);
            obj.layersMap = ic;          
        end

        function [rigidBodySetsPerLayer, numRigids] = buildLayers(obj, mesh3D, cache, h)
            if obj.LayerOrderingType == 1 %random
                elementsArray = 1:size(mesh3D.t,1);
                strainInds = elementsArray(randperm(length(elementsArray)))';
            else %0 : strain rate finite differences
                Fa = mesh3D.B * mesh3D.p;
                Fb = mesh3D.B * cache.oldp; % TODO: use a cached F as these were computed last step
                EDotNorms = mexEdiffNorm3D( Fa, Fb,h );
                [~,strainInds] = sort(EDotNorms,"ascend");
            end

            %incrementally find connectivity for the layers
            %TODO: handle the pinned verts as joints
            forcedElastic= false(size(mesh3D.t,1),1);
            forcedElastic(mesh3D.pinnedTets) = true;

            [rigidBodySetsPerLayer, numRigids] = mexUnionFindByVertice(strainInds, obj.sortedLayers, forcedElastic, mesh3D.t, mesh3D.N);
            rigidBodySetsPerLayer = int32(rigidBodySetsPerLayer);
        end

        function integrate( obj, mesh3D, h, Jc, phi, cInfo, settings, cache, td, animationScripts, frame, energyModel, collisionDectectors, rigidificator)
            assert(~isempty(obj.layers) && ~isempty(obj.iterations));
            assert(obj.boundaryCompliance < 1e-5);
            isOk = ~settings.plotFrameTimeHistogram || obj.runUntilResSmallerThan > 0 || obj.StopAtImprovementPercent > 0;
            assert(isOk);

            sumNonMonitoredTime = 0;%constraint evaluation for residual monitoring is not optimized so we can ignore it in the benchmarks

            if(any(isnan(cache.oldp)) || frame == 1)
                cache.oldp = mesh3D.p;
                obj.pinnedx = mesh3D.p(mesh3D.pinnedDOFs);
            end
            
            ticIntegrateForces = tic;
           
            %this is needed to prevent the rigidificator from complaining
            %about size
            mesh3D.resetForce;
            mesh3D.f = obj.scriptedForceInjection(mesh3D,animationScripts,frame,h);

            % Gravity (z is vertical in matlab)
            mesh3D.f(3:3:end) = mesh3D.f(3:3:end) + mesh3D.mass(3:3:end)*obj.Gravity;
            mesh3D.f = mesh3D.f  + obj.CustomForce; %+ cache.elasticForces

            %build layers
            [rigidBodySetsPerLayer, numRigids] = buildLayers(obj, mesh3D, cache, h);
            isRigidElement = rigidBodySetsPerLayer ~= 0;
            cache.layers = isRigidElement;

            %finding boundaries through intersection, The 2D function works
            %just as well for this
            [isBoundaryVertex,isElasticElement,isRigidVertex,vertexRigidBody] = xpbdLayer2D.findBoundaryVertices(obj.layers, obj.sortedLayers, mesh3D, isRigidElement, rigidBodySetsPerLayer);

            %setting current and prev positions according to the model
            cache.oldOldp = cache.oldp;
            cache.oldp = mesh3D.p;
            xi = mesh3D.p;

            % doing an explicit step
            % mesh3D.v(mesh3D.ElasticDOFs) = mesh3D.v(mesh3D.ElasticDOFs) + h* mesh3D.f(mesh3D.ElasticDOFs)./mesh3D.mass(mesh3D.ElasticDOFs);
            % xi(mesh3D.ElasticDOFs) = xi(mesh3D.ElasticDOFs)  + h * mesh3D.v(mesh3D.ElasticDOFs) ;
            residualVelocity = mesh3D.v;
            residualVelocity(mesh3D.unpinnedDOFs) = residualVelocity(mesh3D.unpinnedDOFs) + h* mesh3D.f(mesh3D.unpinnedDOFs)./mesh3D.mass(mesh3D.unpinnedDOFs);
            residualVelocity = (1-obj.airDrag).*residualVelocity;

            lambdai = zeros(size(mesh3D.t,1),6);%diag and off diagonals of green strain E

            %various settings and setups to monitor perfs
            ticNonMonitoredTime = tic;
            initialResisual = 0;
            cache.initialResidual = initialResisual;
            sumNonMonitoredTime = toc(ticNonMonitoredTime) + sumNonMonitoredTime;

            residualArray = [cache.initialResidual]; 
            contactNormals = cat(1,cInfo.normal);
            pointColliders = cat(1,cInfo.pointCollider);
            contactVertexID = cat(1,cInfo.vertexID);
            contactCompliance = cat(3,cInfo.xpbdContactCompliance);
            if size(contactVertexID,2) > 0
                contactVertexID = contactVertexID(:,1);
            end

            %project
            [xi, sumNonMonitoredTime, residualArray,obj.runUntilResSmallerThan] = xpbdLayer3D.layerSolve(isElasticElement,...
                                                                                                isRigidElement, ...
                                                                                                rigidBodySetsPerLayer, ...
                                                                                                vertexRigidBody, ...
                                                                                                isBoundaryVertex, ...
                                                                                                isRigidVertex, ...
                                                                                                xi, ...
                                                                                                lambdai, ...
                                                                                                sumNonMonitoredTime, ...
                                                                                                h , ...
                                                                                                numRigids, ...
                                                                                                residualArray, ...
                                                                                                obj.boundaryCompliance,...
                                                                                                obj.contactResistance,...
                                                                                                contactCompliance,...
                                                                                                obj.iterations,...
                                                                                                obj.layers,...
                                                                                                obj.runUntilResSmallerThan,...
                                                                                                obj.computeResiduals,...
                                                                                                obj.hangingStop,...
                                                                                                obj.giveUpEnabled,...
                                                                                                obj.giveUpThreshold,...
                                                                                                obj.LastLayerGiveUpEnabled,...
                                                                                                obj.useGravityConstraints,...
                                                                                                cache.oldp,...
                                                                                                cache.oldOldp,...
                                                                                                contactNormals,...
                                                                                                pointColliders ,...
                                                                                                contactVertexID,...
                                                                                                mesh3D.mass,...
                                                                                                mesh3D.t,...
                                                                                                mesh3D.unpinnedDOFs,...
                                                                                                mesh3D.pinnedDOFs,...
                                                                                                mesh3D.DmInv,...
                                                                                                mesh3D.perElementXPBDalphaMatrix,...
                                                                                                mesh3D.perElementXPBDalphaInvMatrix,...
                                                                                                mesh3D.elementColor,...
                                                                                                mesh3D.elementDOFs,...
                                                                                                mesh3D.numColors,...
                                                                                                mesh3D.f,...
                                                                                                obj.layersMap,...
                                                                                                obj.gravityCompliance, ...
                                                                                                obj.Gravity, ...
                                                                                                residualVelocity, ...
                                                                                                mesh3D.perElementXPBDbeta, ...
                                                                                                mesh3D.elV,...
                                                                                                obj.runOnGPU,...
                                                                                                obj.useRigidContacts,...
                                                                                                obj.StopAtImprovementPercent,...
                                                                                                obj.pinnedx);
            
            cache.residualArray = residualArray(:,2:end);

            mesh3D.v = (xi-cache.oldp)./h;     
            mesh3D.p = xi;

            mesh3D.v = xpbdLayer2D.contactVelocityUpdate(contactNormals, contactVertexID,  mesh3D.v, h, obj.Gravity, 3, 1e-2); %The 2D code is compatible here too

            td.integrateForces = toc( ticIntegrateForces ) - sumNonMonitoredTime;
            td.residual = sum(residualArray(:,end));
            td.residualPreSolve = sum(residualArray(:,1));
        end
    end
end