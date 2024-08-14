classdef xpbdResLayer2D < xpbdLayer2D
    % xpbd time integrator. As of now it assumes that you are using solely
    % adaptive meshes

    properties
        LayerPropagation = -1; %-1: all   0:none  1: one
    end
    
    methods(Static)
        function residualv = rotateResidual(rigidIDbyVert, residualv, Rstack) 
            %Rstack is a 3D array of rotation matrices 2by2byNrigids
            rigidVerts = find(rigidIDbyVert ~= 0);
            resvFormatted = formatPositions2D(residualv);
            resvCropped = resvFormatted(rigidVerts,:);
            newResv = pagemtimes(Rstack(:,:,rigidIDbyVert(rigidVerts)),permute(resvCropped,[2,3,1]));
            residualv([rigidVerts*2-1,rigidVerts*2]) = permute(newResv,[3,1,2]);
        end

        function[xi, sumNonMonitoredTime, residualArray, runUntilResSmallerThan]=layerSolve(isElasticElement, ...
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
                                                                                    perElementXPBDBeta, ...
                                                                                    runOnGPU,...
                                                                                    setPercentStop,...
                                                                                    volume,...
                                                                                    useRigidConstacts)
            errorChangeDecreased = false;
            toGPUtic = tic;
            if runOnGPU
                [residualv, isElasticElement,isRigidElement,rigidBodySetsPerLayer,vertexRigidBody,isBoundaryVertex,isRigidVertex, xi, lambdai, h , numRigids, residualArray, boundaryCompliance, contactResistance, contactCompliance, iterations, layers, runUntilResSmallerThan, computeResiduals, hangingStop, giveUpThreshold, oldp, oldOldp, contactNormals, pointColliders, contactVertexID, mass, T, unpinnedDOFs, pinnedDOFs, DmInv, perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix, elementColor, elementDOFs, numColors, f, layersMap, gravityCompliance, gravity, perElementXPBDBeta]=toGPUArray(residualv, isElasticElement, isRigidElement, rigidBodySetsPerLayer, vertexRigidBody, isBoundaryVertex, isRigidVertex, xi, lambdai, h , numRigids, residualArray, boundaryCompliance, contactResistance, contactCompliance, iterations, layers, runUntilResSmallerThan, computeResiduals, hangingStop, giveUpThreshold, oldp, oldOldp, contactNormals, pointColliders, contactVertexID, mass, T, unpinnedDOFs, pinnedDOFs, DmInv,perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix, elementColor, elementDOFs, numColors, f, layersMap, gravityCompliance, gravity, perElementXPBDBeta);
                lambdaGravity = zeros(numel(unpinnedDOFs)/2,1,'gpuArray');
                contactLambdas = zeros(numel(contactVertexID),1,'gpuArray');
            else
                lambdaGravity = zeros(numel(unpinnedDOFs)/2,1);
                contactLambdas = zeros(numel(contactVertexID),1);
            end
            sumNonMonitoredTime = sumNonMonitoredTime + toc(toGPUtic);

            LayerPropagation =-1;

            coder.gpu.nokernel();
            xprev = xi;
            for layer = 1:numel(layers)
                isLastLayer = layer == numel(layers);
                mappedLayerIndex = layersMap(layer);
                elasticElements = find(isElasticElement(:,mappedLayerIndex));
                rigidElements = find(isRigidElement(:,mappedLayerIndex));
                rigidIDbyElement = rigidBodySetsPerLayer(:,mappedLayerIndex);
                rigidIDbyVert = vertexRigidBody(:,mappedLayerIndex);
                isVertBoundary = isBoundaryVertex(:,mappedLayerIndex);
                rigidVerts = find(isRigidVertex(:,mappedLayerIndex));
                elasticVerts = find(~isRigidVertex(:,mappedLayerIndex));
                ElasticDOFs = reshape([elasticVerts*2-1;elasticVerts*2],[],1);

                %build rigids from components
                [R, COM, angle, inertia, rigidMasses, ri, COMDOT, angularv]=AdaptiveMesh.makeRigidsFromConnectivity(numRigids(mappedLayerIndex), xi, residualv, rigidIDbyVert, mass);
                if numel(angle) == 0 && numel(rigidElements)>0
                    disp("failure to build any rigid body");
                end

                %step rigids and released particles
                %release all elastic residual
                xi(ElasticDOFs) = xi(ElasticDOFs) + h*residualv(ElasticDOFs);
                residualv(ElasticDOFs) = 0;
                steppedVelocityR = zeros(2,2,numRigids(mappedLayerIndex));

                if (layer <= LayerPropagation) || (LayerPropagation == -1 && numel(rigidElements)>0)
                    isRigidVert = rigidIDbyVert > 0;
                    xprev = xi;

                    %step rigids
                    [steppedVelocityR, COM, angle] = RigidBody.rigidSymplecticImperativeStep(COM, angle, h, COMDOT, angularv);

                    %update pos wrt stepped rigids
                    xi = rotateAllVerts2D(rigidIDbyVert, ri, steppedVelocityR, COM, xi);
                    
                    rigidInds = find(isRigidVert);
                    rigidDOFs = reshape([rigidInds*2-1, rigidInds*2]',[],1);
                    rigidResidualv = zeros(size(xi,1),1);
                    rigidResidualv(rigidDOFs) = BDF.BDF1(xi(rigidDOFs), xprev(rigidDOFs),h);

                    rotatedResv = rotateAllVerts2D(rigidIDbyVert,formatPositions2D(residualv),steppedVelocityR,zeros(1,2,size(steppedVelocityR,3)),zeros(size(xi)));

                    residualv(rigidDOFs) = rotatedResv(rigidDOFs) - rigidResidualv(rigidDOFs);

                    % Perhaps we should add this to ignore residuals
                    % smaller than numerical accuracy
                    % residualv(abs(residualv)<(h*h)) = 0;
                end

                [xi, lambdai,  lambdaGravity, residualArray,sumNonMonitoredTime,errorChangeDecreased, stopSolve, projectedR, COM, angle, inertia, rigidMasses, runUntilResSmallerThan] = xpbdLayer2D.projectGS(xi, lambdai, lambdaGravity, h, elasticElements, isVertBoundary, residualArray, sumNonMonitoredTime,errorChangeDecreased, contactLambdas,boundaryCompliance,...
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
                                                                                                                                                        angle, ...
                                                                                                                                                        inertia, ...
                                                                                                                                                        rigidMasses, ...
                                                                                                                                                        ri,...
                                                                                                                                                        gravityCompliance, ...
                                                                                                                                                        gravity, ...
                                                                                                                                                        perElementXPBDBeta, ...
                                                                                                                                                        runOnGPU,...
                                                                                                                                                        setPercentStop,...
                                                                                                                                                        volume,...
                                                                                                                                                        useRigidConstacts,...
                                                                                                                                                        norm(residualv));
                if (layer <= LayerPropagation) || (LayerPropagation == -1 && numel(rigidElements)>0)
                    constraintedR = cat(3,projectedR);
                    deltaR = pagemtimes(constraintedR,'none',steppedVelocityR,'transpose');
                    residualv = xpbdResLayer2D.rotateResidual(rigidIDbyVert, residualv, deltaR);
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
    end
    methods
        function obj = xpbdResLayer2D()
            obj@xpbdLayer2D();
            obj.Name = 'xpbdResLayer2D';
        end

        function [mesh2d, cache, td]=integrate( obj, mesh2d, h, Jc, phi, cInfo, settings, cache, td, animationScripts, frame  )
            if nargin < 4
                Jc = zeros( 0, mesh2d.N*2 ); % no constraints
            end          
            if nargin < 5
                phi = [];
            end

            if(isempty(obj.alphaTildeMatrices))
                obj.alphaTildeMatrices = mesh2d.perElementXPBDalphaMatrix./(h*h);
            end
            
            sumNonMonitoredTime = 0;
            ticIntegrateForces = tic;
            
            mesh2d.resetForce;
            mesh2d.f = obj.scriptedForceInjection(mesh2d,animationScripts,frame,h);
            mesh2d.f(2:2:end) = mesh2d.f(2:2:end) + mesh2d.mass(2:2:end)*obj.Gravity;
            mesh2d.f = mesh2d.f  + obj.CustomForce; %+ cache.elasticForces

            sumNonMonitoredTime = 0;

            %incrementally build the connectivity for the layers ------
            ticConnectedComponents = tic;
            [rigidBodySetsPerLayer, numRigids] = buildLayers(obj,mesh2d, cache, h);
            td.connectedComponents = toc(ticConnectedComponents);
            isRigidElement = rigidBodySetsPerLayer ~= 0;
            cache.layers = isRigidElement;

            %finding boundaries through intersection
            [isBoundaryVertex,isElasticElement,isRigidVertex,vertexRigidBody] = xpbdLayer2D.findBoundaryVertices(obj.layers, obj.sortedLayers, mesh2d, isRigidElement, rigidBodySetsPerLayer);

            xi = mesh2d.p;

            % doing an explicit step for velocities only---------------
            residualVelocity = mesh2d.v;
            residualVelocity(mesh2d.unpinnedDOFs) = residualVelocity(mesh2d.unpinnedDOFs) + h* mesh2d.f(mesh2d.unpinnedDOFs)./mesh2d.mass(mesh2d.unpinnedDOFs);


            %various settings and setups to monitor perfs
            cache.xarray = xi;

            lambdai = zeros(size(mesh2d.t,1),3);
            ticNonMonitoredTime = tic;

            initialResisual = 0;
            cache.initialResidual = initialResisual;
            sumNonMonitoredTime = toc(ticNonMonitoredTime) + sumNonMonitoredTime;
            %pinned vertices are hard constraints so we ignore them in the
            %residual
            residualArray = [initialResisual]; 

            
            contactNormals = cat(1,cInfo.normal);
            pointColliders = cat(1,cInfo.pointCollider);
            contactVertexID = cat(1,cInfo.vertexID);

            %project-------------------------
            [xi, sumNonMonitoredTime, residualArray, obj.runUntilResSmallerThan] = xpbdResLayer2D.layerSolve(isElasticElement,...
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
                                                                                                obj.contactCompliance,...
                                                                                                obj.iterations,...
                                                                                                obj.layers,...
                                                                                                obj.runUntilResSmallerThan,...
                                                                                                obj.computeResiduals,...
                                                                                                obj.hangingStop,...
                                                                                                obj.giveUpEnabled,...
                                                                                                obj.giveUpThreshold,...
                                                                                                obj.LastLayerGiveUpEnabled,...
                                                                                                obj.useGravityConstraints,...
                                                                                                mesh2d.prevp,...
                                                                                                mesh2d.prevPrevp,...
                                                                                                contactNormals,...
                                                                                                pointColliders ,...
                                                                                                contactVertexID,...
                                                                                                mesh2d.mass,...
                                                                                                mesh2d.t,...
                                                                                                mesh2d.unpinnedDOFs,...
                                                                                                mesh2d.pinnedDOFs,...
                                                                                                mesh2d.DmInv,...
                                                                                                mesh2d.perElementXPBDalphaMatrix,...
                                                                                                mesh2d.perElementXPBDalphaInvMatrix,...
                                                                                                mesh2d.elementColor,...
                                                                                                mesh2d.elementDOFs,...
                                                                                                mesh2d.numColors,...
                                                                                                mesh2d.f,...
                                                                                                obj.layersMap,...
                                                                                                obj.gravityCompliance, ...
                                                                                                obj.Gravity,...
                                                                                                residualVelocity, ...
                                                                                                mesh2d.perElementXPBDbeta,...
                                                                                                obj.runOnGPU,...
                                                                                                obj.StopAtImprovementPercent,...
                                                                                                mesh2d.elA,...
                                                                                                obj.useRigidConstacts);
            cache.residualArray = xpbdLayer2D.prepResidual(settings,residualArray);
            cache.initialResidual = cache.residualArray(1);

            %update velocities and positions-------------------------
            mesh2d.v = BDF.BDF1(xi,mesh2d.prevp,h);
            % mesh2d.v = BDF.BDF2(xi,mesh2d.prevp,mesh2d.prevPrevp,h);
            % mesh2d.v = BDF.BDF3(xi,mesh2d.prevp,mesh2d.prevPrevp, mesh2d.prevPrevPrevp,h);
            mesh2d.p = xi;

            mesh2d.v = xpbdLayer2D.contactVelocityUpdate(contactNormals, contactVertexID, mesh2d.v, h, obj.Gravity);

            td.residual = sum(residualArray(:,end));
            td.residualPreSolve = sum(residualArray(:,1));
            td.integrateForces = toc( ticIntegrateForces ) - sumNonMonitoredTime;%TODO: remove time from info gathering
        end
    end
end