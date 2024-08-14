classdef xpbdLayer2D < Integrator
    % xpbd time integrator. As of now it assumes that you are using solely
    % adaptive meshes

    %TODO: I should clean this up at some point

    properties
        % none needed
        prevWarmstarts
        regularizator = 1;
        useGamma = 0;
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
        giveUpThreshold = 1e-8;
        contactCompliance = 1e-4;
        contactResistance = 1;
        gravityCompliance = 0.001;
        residualType = 0; %0: constraint residual 1: eq4 residual 2:Ax-b
        LayerOrderingType = 0; %0: strain rate 1: random increment 2: eigs MinvK 3: element order 4:custom
        runUntilResSmallerThan = -1; %if <=0: stops at the given number of iterations. if >0: tries to reach desired residual specified by this variable
        StopAtImprovementPercent = 0;
        hangingStop = 1000; %if runUntil does more than this number of iterations it can be considered as hanging so we stop
        useGravityConstraints = true;
        customOrdering = [];
        alphaTildeMatrices = [];
        useRigidConstacts = true;

        isLastLayer = false;
        runOnGPU = true;

        %for newton residual monitoring
        A;
        b;
        
    end
    properties (Constant = true)%add the energy type here so the parallel loop does not hang as much
        kernelSize = [1,690,1];
    end
    
    methods(Static)
        
        function force = scriptedForceInjection(mesh2d,animationScripts,frame,h)
            force = mesh2d.f;
            for animID = 1:numel(animationScripts)
                force = animationScripts{animID}.scriptForceAnimation(force,frame,h);
            end
        end

        function [lambdai,xi] = stvkVoigtConstraints(numColors, elementDOFs, elementColor, DmInv, alphaTilde, mass, lambdai,xi,elasticElements, betaTilde, oldp, h, volume) %#codegen
            graphColors = elementColor(elasticElements);

            minv = 1./mass;
             %https://www.mathworks.com/help/gpucoder/ref/coder.gpu.nokernel.html
            coder.gpu.nokernel();
            for c = 1 : numColors
                elementsOfColorLogical = graphColors==c;
                elementsOfColor = elasticElements(elementsOfColorLogical);                
                colorDOFs = elementDOFs(elementsOfColor,:);

                xiOfColor = reshape(xi(colorDOFs)',6,1,[]);
                xiOldOfColor = reshape(oldp(colorDOFs)',6,1,[]);
                minvOfColor = reshape(minv(colorDOFs)',6,1,[]);
                lambdaOfColor = reshape(lambdai(elementsOfColor,:)',3,1,[]);
                alphaTildeOfColor = alphaTilde(:,:,elementsOfColor);
                betaTildeOfColor = betaTilde(:,:,elementsOfColor);
                DmInvOfColor = DmInv(:,:,elementsOfColor);
                volumeOfColor = reshape(volume(elementsOfColor),1,1,[]);

                [lambdaOfColor,xiOfColor] = StVKForColor(lambdaOfColor, xiOfColor, xiOldOfColor, alphaTildeOfColor, minvOfColor, DmInvOfColor, betaTildeOfColor, h, volumeOfColor);
                % [lambdaParallel,xiParallel] = StVKForColor_mex(elementsOfColor, lambdaParallel, colorDOFs, xi, xiParallel, alpha, mesh2d.mass, mesh2d.DmInv, h);

                %This overwritting is only ok because the colors are
                %hitting independant dofs
                lambdai(elementsOfColor,:)=reshape(lambdaOfColor,3,[])';
                xi(colorDOFs) = reshape(xiOfColor,6,[])';
            end
        end

        function [lambdaGravity, xi, constraintC] = gravityConstraints(mass, unpinnedDOFs, lambdaGravity, xi, h, prevx, gravityCompliance, gravity) %#codegen

            dofs = formatPositions2D(unpinnedDOFs);
            parallelxi = xi(dofs(:,2));
            parallelxold = prevx(dofs(:,2));
            parallelMass = mass(dofs(:,2));

            coder.gpu.kernelfun();
            [deltax, deltaLambda, constraintC] = gravityConstraint.project(parallelxi, parallelxold, lambdaGravity, parallelMass, gravityCompliance , gravity ,h);
            lambdaGravity = lambdaGravity + deltaLambda;
            xi(dofs(:,2)) = parallelxi + deltax;           
        end


        function [xi, lambdai,  lambdaGravity, residualArray,sumNonMonitoredTime,errorChangeDecreased, stopSolve, R, COM, angle, inertia, rigidMasses, runUntilResSmallerThan] = projectGS(xi, ...
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
                                                                                                                                                        angle, ...
                                                                                                                                                        inertia, ...
                                                                                                                                                        rigidMasses, ...
                                                                                                                                                        ri,...
                                                                                                                                                        gravityCompliance, ...
                                                                                                                                                        gravity, ...
                                                                                                                                                        perElementXPBDBeta, ...
                                                                                                                                                        runOnGPU, ...
                                                                                                                                                        setPercentStop,...
                                                                                                                                                        volume,...
                                                                                                                                                        useRigidConstacts,...
                                                                                                                                                        resvnorm)
    
            
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
                [lambdai,xi] = xpbdLayer2D.stvkVoigtConstraints(numColors, elementDOFs, elementColor, DmInv, perElementXPBDalphaMatrix./(h*h), mass, lambdai,xi,elasticElements, perElementXPBDBeta.*(h*h), oldp, h, volume);

                %rigid-elastic glue distance constraints 
                [lambdaBoundary, xi, R, COM, angle, inertia, rigidMasses] = xpbdLayer2D.rigidElasticBoundaryAccumulationConstraints(boundaryCompliance, mass, rigidIDbyVert, R, COM, angle, inertia, rigidMasses, ri, xi, boundaryVertToSolve, lambdaBoundary, h); 

                %to handle contacts as hard constraints
                [contactLambdas, xi, R, COM, angle, inertia, rigidMasses] = xpbdLayer2D.contactConstraints(contactNormals, pointColliders, contactVertexID, mass, rigidIDbyVert, contactLambdas, xi, h, contactResistance, contactCompliance,R, COM, angle, inertia, rigidMasses, ri, useRigidConstacts);

                %handling pinned constraints
                xi(pinnedDOFs) = oldp(pinnedDOFs);

                ticNonMonitoredTime = tic;

                [endLayer, residualArray, errorChangeDecreased, ~, stopSolve, runUntilResSmallerThan, setPercentStop] = xpbdLayer2D.simpleMonitorError(h, lambdai, xi, oldp, oldOldp, mass, T, f, pinnedDOFs, elementDOFs, DmInv, perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix, residualArray, errorChangeDecreased, computeResiduals, runUntilResSmallerThan, hangingStop, giveUpEnabled, giveUpThreshold, isLastLayer, LastLayerGiveUpEnabled, lambdaGravity, gravityCompliance, gravity, unpinnedDOFs, contactNormals, pointColliders, contactVertexID, contactLambdas, contactResistance, contactCompliance, ri, perElementXPBDBeta, setPercentStop, volume, resvnorm );

                sumNonMonitoredTime =  toc(ticNonMonitoredTime) + sumNonMonitoredTime;
                if endLayer || stopSolve
                    break;
                end
            end

            %because no rigid body share dofs, we can vectorize
            xi = rotateAllVerts2D(rigidIDbyVert, ri, R, COM, xi);
        end

        function residual = computePositionLevelResidual(M,x,xprev,xprevprev,alpha,gradC,C,h)
            %residual inspired by XPBD paper equation 4
            h2 = h*h;
            midx = (x - 2*xprev + xprevprev)./h2;
            Lhs = M*midx;

            rhs = -gradC*(alpha\C);

            residual = norm(Lhs - rhs);
        end

        function residual = computeResidual(constraintC,gradC,Minv,alphaTilde,lambda, deltaLambda)
            residual = (gradC'*Minv*gradC + alphaTilde)*deltaLambda + constraintC + alphaTilde*lambda;  
        end

        function residual = computeResidualh(constraintC,alphaTilde,lambda)
            residual = constraintC + alphaTilde*lambda;  %ignore damping... for now I set it to 0
        end


        function [lambdai,xi] = stvkConstraints(obj, mesh2d, lambdai,xi,elasticElements, cache, h)
            elementDOFs = mesh2d.elementDOFs;
            graphColors = mesh2d.elementColor(elasticElements);

            coder.gpu.nokernel();
            for c = 1 : mesh2d.numColors
                coder.gpu.kernelfun();

                elementsOfColorLogical = graphColors==c;
                elementsOfColor = elasticElements(elementsOfColorLogical);
                lambdaParallel = lambdai(elementsOfColor);
                colorDOFs = elementDOFs(elementsOfColor,:);
                xiParallel = xi(elementDOFs(elementsOfColor,:));
                if size(xiParallel,2)<6
                    xiParallel = xiParallel';
                end
                
                for elasticInd = 1: numel(elementsOfColor)
                    element = elementsOfColor(elasticInd);
                    tetDofsIdsVector = colorDOFs(elasticInd,:)';
    
                    %matrix alpha version:mesh2d.perElementXPBDalphaMatrix(:,:,element)
                    %[deltaLambdai, deltaxi, forces, alphaTilde, constraintC, gradC, Minv] = stvkConstraint.project(mesh2d.prevp(tetDofsIdsVector), xi(tetDofsIdsVector), mesh2d.DmInv(:,:,element), obj.elasticCompliance, obj.elasticBeta, lambdai(element), mesh2d.elMu(element), mesh2d.elLambda(element), mesh2d.mass(tetDofsIdsVector), mesh2d.elA(element), h);
                    [deltaxi, Minv, gradC, deltaLambdai, alphaTilde, constraintC] = mexSTVKConstraint( h, xi(tetDofsIdsVector), mesh2d.prevp(tetDofsIdsVector), mesh2d.DmInv(:,:,element), mesh2d.elMu(element), mesh2d.elLambda(element), mesh2d.elA(element), obj.elasticCompliance, obj.elasticBeta, lambdaParallel(elasticInd), mesh2d.mass(tetDofsIdsVector) );
                    
                    %I expect that only computing the residual on all the elastic
                    %constraints is enough.... we don't really care about the rest
                    %as they are transitionary
                    lambdaParallel(elasticInd) = lambdaParallel(elasticInd) + deltaLambdai;
                    xiParallel(elasticInd,:) = xiParallel(elasticInd,:) + deltaxi';
                end
                lambdai(elementsOfColor)=lambdaParallel;
                xi(colorDOFs) = xiParallel;
            end
        end


        function [lambdai,xi] = neoHookeanConstraints(obj, mesh2d, lambdai,xi,elasticElements, h)
            elementDOFs = mesh2d.elementDOFs;
            graphColors = mesh2d.elementColor(elasticElements);

            coder.gpu.nokernel();
            for c = 1 : mesh2d.numColors
                coder.gpu.kernelfun();
                elementsOfColorLogical = graphColors==c;
                elementsOfColor = elasticElements(elementsOfColorLogical);
                lambdaParallelD = lambdai(elementsOfColor,1);
                lambdaParallelH = lambdai(elementsOfColor,2);
                colorDOFs = elementDOFs(elementsOfColor,:);
                xiParallel = xi(elementDOFs(elementsOfColor,:));
                if size(xiParallel,2)<6
                    xiParallel = xiParallel';
                end

                
                for elasticInd = 1: numel(elementsOfColor)%TODO: back to parallel?
                    element = elementsOfColor(elasticInd);
                    tetDofsIdsVector = colorDOFs(elasticInd,:)';
    
                    [deltaxi, Minv, gradC, deltaLambdai, alphaTildeD, constraintCD] = mexNeoHookeanDConstraint(h, xi(tetDofsIdsVector), mesh2d.DmInv(:,:,element), mesh2d.elMu(element), mesh2d.elLambda(element), mesh2d.elA(element), lambdaParallelD(elasticInd), mesh2d.mass(tetDofsIdsVector));
                    lambdaParallelD(elasticInd) = lambdaParallelD(elasticInd) + deltaLambdai;
                    xiParallel(elasticInd,:) = xiParallel(elasticInd,:) + deltaxi';

                    [deltaxi, Minv, gradC, deltaLambdai, alphaTildeH, constraintCH] = mexNeoHookeanHConstraint(h, xi(tetDofsIdsVector), mesh2d.DmInv(:,:,element), mesh2d.elMu(element), mesh2d.elLambda(element), mesh2d.elA(element), lambdaParallelH(elasticInd), mesh2d.mass(tetDofsIdsVector), mesh2d.elGamma(element));
                    lambdaParallelH(elasticInd) = lambdaParallelH(elasticInd) + deltaLambdai;
                    xiParallel(elasticInd,:) = xiParallel(elasticInd,:) + deltaxi';
                end
                lambdai(elementsOfColor,1)=lambdaParallelD;
                lambdai(elementsOfColor,2)=lambdaParallelH;
                xi(colorDOFs) = xiParallel;
            end
        end

        function [lambdaBoundary, xi, R, COM, angle, inertia, rigidMasses] = rigidElasticBoundaryAccumulationConstraints(boundaryCompliance, mass, rigidIDbyVert, R, COM, angle, inertia, rigidMasses, ri, xi, boundaryVertToSolve, lambdaBoundary, h)%#codegen
            if(numel(lambdaBoundary)>0)
                bodyIDs = rigidIDbyVert(boundaryVertToSolve);
                boundaryElasticDOFs = reshape([boundaryVertToSolve*2-1,boundaryVertToSolve*2]',1,2,[]);
                elasticxi = xi(boundaryElasticDOFs);
                massInvElastic = reshape(1./mass(boundaryVertToSolve*2),1,1,[]);
    
                riCropped = permute(ri(boundaryVertToSolve,:),[3,2,1]);
                rigidxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,bodyIDs), COM(:,:,bodyIDs));
                rigidInertia = reshape(inertia(bodyIDs),1,1,[]);
                rigidAngle = angle(bodyIDs);
    
                r = rigidxi-COM(:,:,bodyIDs);
                boundaryVertsCount = accumarray(bodyIDs,ones(size(bodyIDs)));
                splitUpdateByVerts = reshape(boundaryVertsCount(bodyIDs),1,1,[]);
    
                [deltaLambda, ~, deltaCOM, deltaOmega]=rigidBoundaryConstraint2D.projectVectorized(h, elasticxi, rigidxi, r, massInvElastic, reshape(lambdaBoundary,1,1,[]), boundaryCompliance, reshape(rigidMasses(bodyIDs),1,1,[]), pagetranspose(COM(:,:,bodyIDs)), rigidAngle, rigidInertia, splitUpdateByVerts);
    
                lambdaBoundary = lambdaBoundary + reshape(deltaLambda,[],1,1);
                numRigids = size(R,3);

                deltaCOMxAccumulator = reshape(accumarray(bodyIDs,reshape(deltaCOM(:,1,:),[],1,1),[numRigids,1]),1,1,[]);
                deltaCOMyAccumulator = reshape(accumarray(bodyIDs,reshape(deltaCOM(:,2,:),[],1,1),[numRigids,1]),1,1,[]);
                deltaCOMAccumulator = cat(2,deltaCOMxAccumulator,deltaCOMyAccumulator);
                deltaOmegaAccumulator = accumarray(bodyIDs,reshape(deltaOmega,[],1,1),[numRigids,1]);

                COM = COM + deltaCOMAccumulator;
                angle(:) = angle(:) + deltaOmegaAccumulator;
                R = angle2rotm2D(angle);

                riCropped = permute(ri(boundaryVertToSolve,:),[3,2,1]);
                rigidxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,bodyIDs), COM(:,:,bodyIDs));
                xi(boundaryElasticDOFs) = rigidxi; %constraint is exactly enforced to preserve the rigid boundary
            end
        end

        function [lambdaBoundary, xi, R, COM, angle, inertia, rigidMasses] = rigidElasticBoundaryConstraints(boundaryCompliance, mass, rigidIDbyVert, R, COM, angle, inertia, rigidMasses, ri, xi, boundaryVertToSolve, lambdaBoundary, h)%#codegen
            %reorder this to paralellize the rigids bodies by doing the
            %parfor on rigids and sequentially solve constraints within
            %each of them
            numRigids = numel(angle);
            if numRigids == 0
                return
            end

            bodyIDs = rigidIDbyVert(boundaryVertToSolve);
            dofs = [boundaryVertToSolve*2-1,boundaryVertToSolve*2];
           
            listRigidIDs = 1:numRigids;
            xiParallel = cell(numRigids,1);
            lambdaParallel = cell(numRigids,1);
            boundaryVerts = cell(numRigids,1);
            boundaryVerticesOfBodyLogical = (bodyIDs == listRigidIDs)';

            % coder.gpu.kernel(xpbdLayer2D.kernelSize);
            coder.gpu.kernelfun();
            for b = 1:numRigids %TODO: back to parallel
                boundaryVerts{b} = boundaryVertToSolve(boundaryVerticesOfBodyLogical(b,:));
                lambdaParallel{b} = lambdaBoundary(boundaryVerticesOfBodyLogical(b,:));
                xiParallel{b} = xi(dofs(boundaryVerticesOfBodyLogical(b,:),:));
            end

            % coder.gpu.kernel(xpbdLayer2D.kernelSize);
            coder.gpu.kernelfun();
            parfor b = 1:numRigids %TODO: back to parallel
                coder.gpu.nokernel();
                for i = 1 : numel(boundaryVerts{b})
                    vertID = boundaryVerts{b}(i);
                    % vertexDOFs = [vertID*2-1,vertID*2];

                    rigidxi = RigidBody.relativeVertexPosition(ri(vertID,:), R(:,:,b), COM(:,b)');
                    r = rigidxi-COM(:,b)';

                    % [deltaLambda, deltaxElastic, rigidCOM, rigidAngle] = mexRigidBoundaryConstraint(h,xiParallel{b}(i,:)',rigidxi,r, mass(vertID*2), lambdaParallel{b}(i), boundaryCompliance, rigidMasses(b),COM(:,b)', angle(b), inertia(b), 1);
                    [deltaLambda, deltaxElastic, rigidCOM, rigidAngle]=rigidBoundaryConstraint2D.project(h, xiParallel{b}(i,:), rigidxi, r, mass(vertID*2), lambdaParallel{b}(i), boundaryCompliance, rigidMasses(b), COM(:,b)', angle(b), inertia(b), 1);

                    angle(b) = rigidAngle;
                    COM(:,b) = rigidCOM';
                    R(:,:,b) = angle2rotm2D(angle(b));

                    lambdaParallel{b}(i) = lambdaParallel{b}(i) + deltaLambda;
                    xiParallel{b}(i,:)  = xiParallel{b}(i,:) + deltaxElastic;
                end
            end

            for b = 1:numRigids
                xi(dofs(boundaryVerticesOfBodyLogical(b,:),:)) = xiParallel{b};
                lambdaBoundary(boundaryVerticesOfBodyLogical(b,:)) = lambdaParallel{b};
            end
        end


        function [contactLambdas, xi, R, COM, angle, inertia, rigidMasses] = contactConstraints(contactNormals, pointColliders, vertexIDs, mass, rigidIDbyVert, contactLambdas, xi, h, contactResistance, contactCompliance, R, COM, angle, inertia, rigidMasses, ri, useRigidConstacts)%#codegen
            %I probably want to paralellize these... should split between
            %rigid and elastic, then solve the rigid ones sequentially when
            %they belong to the same rigid body
            minv = 1./mass;

            isElasticContact = rigidIDbyVert(vertexIDs) == 0;
            elasticContactVerts = vertexIDs(isElasticContact);
            elasticContactPoints = reshape(pointColliders(isElasticContact,:)',1,2,[]);
            elasticContactNormals = reshape(contactNormals(isElasticContact, :)',1,2,[]);
            elasticContactDOFs = reshape([elasticContactVerts*2-1, elasticContactVerts*2]',1,2,[]); 
            formattedxi = reshape(xi(elasticContactDOFs),1,2,[]);
            elasticminv = reshape(minv(elasticContactDOFs) ,1,2,[]);

            [dx, deltalambda,~]=penaltyConstraint2D.projectElastic(formattedxi, reshape(contactLambdas(isElasticContact),1,1,[]), elasticminv, elasticContactPoints, elasticContactNormals, contactResistance, contactCompliance, h);
            xi(elasticContactDOFs) = formattedxi + dx;
            contactLambdas(isElasticContact) = contactLambdas(isElasticContact) + deltalambda(:);       

            isRigidContact = ~isElasticContact;
            rigidContactVerts = vertexIDs(isRigidContact);
            if numel(rigidContactVerts) <1 || ~useRigidConstacts
                return;
            end
            
            bodyIDs = rigidIDbyVert(rigidContactVerts);
            rigidParticleDOFs = reshape([rigidContactVerts*2-1,rigidContactVerts*2]',1,2,[]);
            rigidParticlexi = reshape(xi(rigidParticleDOFs),1,2,[]);

            riCropped = permute(ri(rigidContactVerts,:),[3,2,1]);
            rigidxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,bodyIDs), COM(:,:,bodyIDs));
            rigidInertia = reshape(inertia(bodyIDs),1,1,[]);
            rigidAngle = angle(bodyIDs);

            r = rigidxi-COM(:,:,bodyIDs);
            boundaryVertsCount = accumarray(bodyIDs,ones(size(bodyIDs)),[size(R,3),1]);
            splitUpdateByVerts = reshape(boundaryVertsCount(bodyIDs),1,1,[]);
            rigidContactPoint = reshape(pointColliders(isRigidContact,:)',1,2,[]);
            rigidContactNormal = reshape(contactNormals(isRigidContact,:)',1,2,[]);
            rigidParticleminv = reshape(minv(rigidParticleDOFs) ,1,2,[]);

            [deltaLambda, deltaCOM, deltaOmega, ~]=penaltyConstraint2D.projectRigid(h, rigidParticlexi, rigidContactPoint, rigidContactNormal, contactResistance, r, reshape(contactLambdas(isRigidContact),1,1,[]), contactCompliance,  reshape(rigidMasses(bodyIDs),1,1,[]), rigidParticleminv(:,1,:), rigidInertia,splitUpdateByVerts);

            contactLambdas(isRigidContact) = contactLambdas(isRigidContact) + reshape(deltaLambda,[],1,1);

            deltaCOMxAccumulator = accumarray(bodyIDs,reshape(deltaCOM(:,1,:),[],1,1),[size(R,3),1]);
            deltaCOMyAccumulator = accumarray(bodyIDs,reshape(deltaCOM(:,2,:),[],1,1),[size(R,3),1]);
            deltaOmegaAccumulator = accumarray(bodyIDs,reshape(deltaOmega,[],1,1),[size(R,3),1]);

            COM(1,1,:) = reshape(COM(1,1,:),[],1,1) + deltaCOMxAccumulator;
            COM(1,2,:) = reshape(COM(1,2,:),[],1,1) + deltaCOMyAccumulator;
            angle(:) = angle(:) + deltaOmegaAccumulator;
            R = angle2rotm2D(angle);

            riCropped = permute(ri(rigidContactVerts,:),[3,2,1]);
            rigidxi = RigidBody.pageRelativeVertexPosition(riCropped,R(:,:,bodyIDs), COM(:,:,bodyIDs));
            xi(rigidParticleDOFs) = rigidxi; %constraint is exactly enforced to preserve the rigid boundary
        end


        function [residualVector] = simpleEvaluateAllErrorEQ4(xi, oldp, oldOldp, mass, T, f, pinnedDOFs, h, lambdai, elementDOFs, DmInv, perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix)%#codegen
            gradU = zeros(size(xi));
            h2 = h*h;
            acceleration = (xi - 2*oldp + oldOldp)./h2;
            Lhs = mass.*acceleration;
            alphaTilde = perElementXPBDalphaMatrix./h2;
            minv = 1./mass;

            coder.gpu.kernelfun();              
            %here everything is considered the same color so we can
            %evaluate in parallel
            xiOfColor = reshape(xi(elementDOFs)',6,1,[]);
            minvOfColor = reshape(minv(elementDOFs)',6,1,[]);
            lambdaOfColor = reshape(lambdai(1:size(T,1),:),3,1,[]);
            [~, gradC, ~, constraintC] = matStvkVoigtConstraint( xiOfColor, xiOldOfColor, DmInv, alphaTilde, lambdaOfColor, minvOfColor, betaTilde,h,1 );

            alphaC = pagemtimes(perElementXPBDalphaInvMatrix,'none',constraintC,'none');
            gradUContributions = pagemtimes(gradC,'transpose',alphaC, 'none');

            coder.gpu.nokernel();
            for elasticInd = 1: size(T,1)%To make sure that the error is fair, we test all elements, but note that the error of rigid elements should be static for a given layer
                gradU(elementDOFs(elasticInd,:)) = gradU(elementDOFs(elasticInd,:)) + gradUContributions(:,:,elasticInd);
            end

            %contact force = lambda *n ./ h^2
            forces = f -gradU;
            forces(pinnedDOFs) = 0;
            residualVector = norm(Lhs-forces,2);
        end

        function [residualVector] = simpleEvaluateAllConstraintResiduals(xi, oldp, oldOldp, mass, T, f, pinnedDOFs, h, lambdai, lambdaGravity, gravityCompliance, gravity, unpinnedDOFs, elementDOFs, DmInv, perElementXPBDalphaMatrix, contactNormals, pointColliders, contactVertexIDs, contactLambdas, contactResistance, contactCompliance, ri, perElementXPBDBeta, volume)
            h2 = h*h;
            alphaTilde = perElementXPBDalphaMatrix./h2;
            betaTilde = perElementXPBDBeta .*h2;
            minv = 1./mass;

            coder.gpu.kernelfun();              
            %here everything is considered the same color so we can
            %evaluate in parallel
            xiOfColor = reshape(xi(elementDOFs)',6,1,[]);
            xiOldOfColor = reshape(oldp(elementDOFs)',6,1,[]);
            minvOfColor = reshape(minv(elementDOFs)',6,1,[]);
            lambdaOfColor = reshape(lambdai(1:size(T,1),:),3,1,[]);
            volumeOfColor = reshape(volume(1:size(T,1)),1,1,[]);
            % volumeOfColor = ones(1,1,size(T,1));
            [~, ~, ~, stvkC] = matStvkVoigtConstraint( xiOfColor, xiOldOfColor, DmInv, alphaTilde, lambdaOfColor, minvOfColor, betaTilde, h, volumeOfColor);
            correctionTerm = abs(pagemtimes(alphaTilde,lambdaOfColor));
            absStvk = abs(stvkC);
            constraintResiduals = absStvk - min(absStvk,correctionTerm);

            %gravity
            % [~, ~, constraintCg] = xpbdLayer2D.gravityConstraints(mass, unpinnedDOFs, lambdaGravity, xi, h, oldp, gravityCompliance, gravity);
            % gravityCorrection = abs((gravityCompliance/h2).*lambdaGravity);
            % absConstraintCg = abs(constraintCg);
            % constraintResidualsGravity = absConstraintCg - min(absConstraintCg,gravityCorrection);


            %here analyses the constraint residual as purely elastic
            %contacts
            contactPoints = reshape(pointColliders',1,2,[]);
            contactNormalsReshaped = reshape(contactNormals',1,2,[]);
            contactDOFs = reshape([contactVertexIDs*2-1, contactVertexIDs*2]',1,2,[]); 
            formattedxi = reshape(xi(contactDOFs),1,2,[]);
            elasticminv = reshape(minv(contactDOFs) ,1,2,[]);

            [~, ~,contactC]=penaltyConstraint2D.projectElastic(formattedxi, reshape(contactLambdas,1,1,[]), elasticminv, contactPoints, contactNormalsReshaped, contactResistance, contactCompliance, h);
            absContact = abs(contactC(:));
            contactCorrection = abs((contactCompliance/h2).*contactLambdas(:));
            constraintResidualsContacts = absContact - min(absContact,contactCorrection);

            residualVector = norm([constraintResiduals(:);constraintResidualsContacts(:)],2); %constraintResidualsGravity;
        end

        function [endLayer, residualArray, errorChangeDecreased, newResidual, stopSolve, runUntilResSmallerThan, setPercentStop] = simpleMonitorError(h, lambdai, xi, oldp, oldOldp, mass, T, f, pinnedDOFs, elementDOFs, DmInv,perElementXPBDalphaMatrix, perElementXPBDalphaInvMatrix, residualArray, errorChangeDecreased, computeResiduals, runUntilResSmallerThan, hangingStop, giveUpEnabled, giveUpThreshold, isLastLayer, LastLayerGiveUpEnabled, lambdaGravity, gravityCompliance, gravity, unpinnedDOFs, contactNormals, pointColliders, contactVertexIDs, contactLambdas, contactResistance, contactCompliance, ri, perElementXPBDBeta, setPercentStop, volume, resvnorm)
            endLayer = false; %move on to next layer
            stopSolve = false; %complete stop of solve
            newResidual = inf;
            if computeResiduals
                newRes = xpbdLayer2D.simpleEvaluateAllConstraintResiduals(xi, oldp, oldOldp, mass, T, f, pinnedDOFs, h, lambdai, lambdaGravity, gravityCompliance, gravity, unpinnedDOFs, elementDOFs, DmInv, perElementXPBDalphaMatrix, contactNormals, pointColliders, contactVertexIDs, contactLambdas, contactResistance, contactCompliance, ri, perElementXPBDBeta, volume);
                residualArray = [residualArray,newRes]; 

                newResidual = norm(newRes);
                lastResidual = norm(residualArray(:,end-1));
                if setPercentStop > 0
                    runUntilResSmallerThan = (1-setPercentStop/100.0)*newResidual;
                    setPercentStop = 0;
                end

                errorChange = newResidual-lastResidual;
                if errorChange < 0 %made this mechanism because sometimes the residual increases at first even in an fully elastic setup... because of this it would make the layer solver give up on layers prematurely
                    errorChangeDecreased = true;
                end

                if giveUpEnabled   && errorChange > giveUpThreshold && errorChangeDecreased && (~isLastLayer || LastLayerGiveUpEnabled) %&& numel(rigidElements) > 0 && it > 1
                    endLayer = true;
                end

                stopHanging = size(residualArray,2) > hangingStop &&  hangingStop >= 0;
                if (runUntilResSmallerThan >= 0 && (runUntilResSmallerThan >= newResidual) && (isLastLayer|| resvnorm < h*h)) || stopHanging
                    stopSolve = true;
                end
            end
        end

        function [isBoundaryVertex,isElasticElement,isRigidVertex,vertexRigidBody] = findBoundaryVertices(layers, sortedLayers, mesh2d, isRigidElement, rigidBodyElementSetsPerLayer)
            isRigidVertex = false(mesh2d.N,numel(layers));

            isElasticElement = ~isRigidElement;
            isVertexOfElasticElement = false(size(isRigidVertex));

            vertexRigidBody = int32(zeros(mesh2d.N,numel(layers)));
            
            for layer = 1:numel(sortedLayers)
                isRigidLayer =  isRigidElement(:,layer);
                isElasticLayer = isElasticElement(:,layer);
                for vertOfT = 1:size(mesh2d.t,2)
                    isRigidVertex(mesh2d.t(isRigidLayer,vertOfT),layer) = true;
                    isVertexOfElasticElement(mesh2d.t(isElasticLayer,vertOfT),layer) = true;
                    vertexRigidBody(mesh2d.t(isRigidLayer,vertOfT),layer) = rigidBodyElementSetsPerLayer(isRigidLayer,layer);
                end
            end
            isBoundaryVertex = isVertexOfElasticElement & isRigidVertex;
        end

        function residualArray = prepResidual(settings, residualArray)
            if settings.skipPreProjectionResidual
                residualArray = residualArray(:,2:end);
            end
        end

        function v = contactVelocityUpdate(contactNormals, contactVertexIDs, v, h, gravity, dim, restitutionCoeff)
            if nargin < 7
                restitutionCoeff = 1e-5;
            end
            if (nargin < 6)
                dim = 2;
            end
            %Here I assume that we can simply consider the contacts to be
            %ultimately solved as elastic. Note that this would potentially
            %create artifacts for purely rigid scenes (could be fixed with a rigid velocity update instead when pure rigidity detected).
            numContacts = size(contactNormals,1);
            if numContacts == 0
                return;
            end
            e = ones(numContacts,1)*restitutionCoeff;
            dofs = repmat(contactVertexIDs*dim,1,dim) - (dim-1:-1:0);
            
            vn = reshape(v(dofs),[],dim);
            vnNorm = vecnorm(vn,2,2);

            dofsUnderRestitutionThreshold = vnNorm <= 2*gravity*h; %taken from the pbdBodies paper
            e(dofsUnderRestitutionThreshold) = 0;

            vTilden = dot(contactNormals,vn,2);
            dv = contactNormals.*(-vn + min(-e.*vTilden,0));
            v(dofs) = vn + dv;
        end

        
        function[xi, sumNonMonitoredTime, residualArray, runUntilResSmallerThan]=layerSolve(isElasticElement,isRigidElement,rigidBodySetsPerLayer,vertexRigidBody,isBoundaryVertex,isRigidVertex, xi, lambdai, sumNonMonitoredTime, h , numRigids, residualArray,boundaryCompliance,...
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
                                                                                                                                                        gravity,...
                                                                                                                                                        perElementXPBDBeta, ...
                                                                                                                                                        runOnGPU, ...
                                                                                                                                                        StopAtImprovementPercent,...
                                                                                                                                                        volume,...
                                                                                                                                                        useRigidConstacts)
            
            errorChangeDecreased = false;

            toGPUtic = tic;
            if runOnGPU
                [isElasticElement,isRigidElement,rigidBodySetsPerLayer,vertexRigidBody,isBoundaryVertex,isRigidVertex, xi, lambdai, sumNonMonitoredTime, h , numRigids, residualArray,boundaryCompliance,contactResistance,contactCompliance,iterations,layers, runUntilResSmallerThan,computeResiduals,hangingStop,giveUpEnabled,giveUpThreshold,LastLayerGiveUpEnabled,useGravityConstraints,oldp,oldOldp,contactNormals,pointColliders,contactVertexID,mass,T,unpinnedDOFs,pinnedDOFs,DmInv,perElementXPBDalphaMatrix,perElementXPBDalphaInvMatrix,elementColor,elementDOFs,numColors,f,layersMap,gravityCompliance, gravity,perElementXPBDBeta]=toGPUArray(isElasticElement,isRigidElement,rigidBodySetsPerLayer,vertexRigidBody,isBoundaryVertex,isRigidVertex, xi, lambdai, sumNonMonitoredTime, h , numRigids, residualArray,boundaryCompliance,contactResistance,contactCompliance,iterations,layers, runUntilResSmallerThan,computeResiduals,hangingStop,giveUpEnabled,giveUpThreshold,LastLayerGiveUpEnabled,useGravityConstraints,oldp,oldOldp,contactNormals,pointColliders,contactVertexID,mass,T,unpinnedDOFs,pinnedDOFs,DmInv,perElementXPBDalphaMatrix,perElementXPBDalphaInvMatrix,elementColor,elementDOFs,numColors,f,layersMap,gravityCompliance, gravity,perElementXPBDBeta);
                lambdaGravity = zeros(numel(unpinnedDOFs)/2,1,'gpuArray');
                contactLambdas = zeros(numel(contactVertexID),1,'gpuArray');
            else
                lambdaGravity = zeros(numel(unpinnedDOFs)/2,1);
                contactLambdas = zeros(numel(contactVertexID),1);
            end
            sumNonMonitoredTime = sumNonMonitoredTime + toc(toGPUtic);

            coder.gpu.nokernel();
            prevp = xi;
            for layer = 1:numel(layers)
                if layer == 1
                    setPercentStop = StopAtImprovementPercent;%convoluted way to get my runUntil variable set to a percentage according to the first residual.
                else 
                    setPercentStop = 0;
                end
                isLastLayer = layer == numel(layers);
                mappedLayerIndex = layersMap(layer);
                elasticElements = find(isElasticElement(:,mappedLayerIndex));
                rigidElements = find(isRigidElement(:,mappedLayerIndex));
                rigidIDbyElement = rigidBodySetsPerLayer(:,mappedLayerIndex);
                rigidIDbyVert = vertexRigidBody(:,mappedLayerIndex);
                elasticVerts = find(rigidIDbyVert~=0);
                elasticDOFs = reshape([elasticVerts*2-1,elasticVerts*2]',[],1);
                isVertBoundary = isBoundaryVertex(:,mappedLayerIndex);
                rigidVerts = find(isRigidVertex(:,mappedLayerIndex));

                %build rigids from components
                vForRigids = (xi-prevp)./h;
                isComponentPinned = false(numRigids(mappedLayerIndex), 1); %TODO: make this compatible (easy, but could take some time to implement).
                [R, COM, angle, inertia, rigidMasses, ri, COMDOT, angularv]=AdaptiveMesh.makeRigidsFromConnectivity(numRigids(mappedLayerIndex), xi, vForRigids, rigidIDbyVert, mass);
                if numel(angle) == 0 && numel(rigidElements)>0
                    disp("failure to build any rigid body");
                end

                [xi, lambdai, lambdaGravity, residualArray,sumNonMonitoredTime,errorChangeDecreased, stopSolve,~, ~, ~, ~, ~, runUntilResSmallerThan] = xpbdLayer2D.projectGS(xi, lambdai, lambdaGravity, h, elasticElements, isVertBoundary, residualArray, sumNonMonitoredTime,errorChangeDecreased, contactLambdas,boundaryCompliance,...
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
                                                                                                                                                        R, ...
                                                                                                                                                        COM, ...
                                                                                                                                                        angle, ...
                                                                                                                                                        inertia, ...
                                                                                                                                                        rigidMasses, ...
                                                                                                                                                        ri,...
                                                                                                                                                        gravityCompliance, ...
                                                                                                                                                        gravity,...
                                                                                                                                                        perElementXPBDBeta,...
                                                                                                                                                        runOnGPU,...
                                                                                                                                                        setPercentStop,...
                                                                                                                                                        volume,...
                                                                                                                                                        useRigidConstacts,...
                                                                                                                                                        inf);

                if stopSolve
                    break;
                end
                prevp = xi;
            end
            
            if runOnGPU
                fromGPUtic = tic;
                [xi, sumNonMonitoredTime, residualArray, runUntilResSmallerThan]=toCPUarray(xi, sumNonMonitoredTime, residualArray, runUntilResSmallerThan);
                sumNonMonitoredTime = sumNonMonitoredTime + toc(fromGPUtic);
            end
        end
    end
    methods
        function obj = xpbdLayer2D()
            obj@Integrator();
            obj.Name = 'xpbdLayer2D';
        end

        function setLayers(obj, layersArray)
            obj.layers = layersArray;
            [obj.sortedLayers, ia, ic] = unique(layersArray);
            obj.layersMap = ic;          
        end

        function [rigidBodySetsPerLayer, numRigids] = buildLayers(obj, mesh2d, cache, h)
            if obj.LayerOrderingType == 1 %random
                elementsArray = 1:size(mesh2d.t,1);
                strainInds = elementsArray(randperm(length(elementsArray)))';
            elseif obj.LayerOrderingType == 2 %eigs
                F = mesh2d.B * mesh2d.p;
                [ii,jj,Cvals,grad,psi]=mexComputeSTVKGradHess2D(F,mesh2d.elA,mesh2d.elMu, mesh2d.elLambda);
                C = sparse(ii,jj,Cvals);
                K = mesh2d.B'*C*mesh2d.B;
                [V,D]=eig(diag(1./mesh2d.mass)*K);
                eigenVals = diag(D);
                freqSum = sum(eigenVals(mesh2d.t),2);
                [sortedStrains,strainInds] = sort(freqSum,'descend');
            elseif obj.LayerOrderingType == 3 %element order (stupid.. don't use other than for debugging purpose)
                strainInds = [1:size(mesh2d.t,1)]';
            elseif obj.LayerOrderingType == 4
                strainInds = obj.customOrdering;
            elseif obj.LayerOrderingType == 5 %strain rate derivative
                F = mesh2d.B * mesh2d.p;
                Fdot = mesh2d.B * mesh2d.v;
                EDotNorms = mexEdotNorm2D( F, Fdot );
                [sortedStrains,strainInds] = sort(EDotNorms,"ascend");
            else %0 : strain rate finite differences
                Fa = mesh2d.B * mesh2d.p;
                Fb = mesh2d.B * mesh2d.prevPrevp; % would be better to cache previous iterations
                EDotNorms = mexEdiffNorm2D( Fa, Fb,h );
                [sortedStrains,strainInds] = sort(EDotNorms,"ascend");
            end

            %incrementally find connectivity for the layers
            %TODO: handle the pinned verts
            forcedElastic= false(size(mesh2d.t,1),1);
            forcedElastic(mesh2d.pinnedTris) = true;

            % [rigidBodySetsPerLayer, numRigids] = mexUnionFindBodies(strainInds, obj.sortedLayers, mesh2d.AdjagencyMatrix, forcedElastic);
            [rigidBodySetsPerLayer, numRigids] = mexUnionFindByVertice(strainInds, obj.sortedLayers, forcedElastic, mesh2d.t, mesh2d.N);
            rigidBodySetsPerLayer = int32(rigidBodySetsPerLayer);
        end

        

        function [mesh2d, cache, td]=integrate( obj, mesh2d, h, Jc, phi, cInfo, settings, cache, td, animationScripts, frame )
            if nargin < 4
                Jc = zeros( 0, mesh2d.N*2 ); % no constraints
            end          
            if nargin < 5
                phi = [];
            end
            assert(~isempty(obj.layers) && ~isempty(obj.iterations));
            assert(obj.boundaryCompliance < 1e-5); %big boundary compliance will cause issues
            isOk = ~settings.plotFrameTimeHistogram || obj.runUntilResSmallerThan > 0 || obj.StopAtImprovementPercent > 0;
            assert(isOk)
            sumNonMonitoredTime = 0;

            if(isempty(obj.alphaTildeMatrices))
                obj.alphaTildeMatrices = mesh2d.perElementXPBDalphaMatrix./(h*h);
            end

            ticIntegrateForces = tic;
            
            mesh2d.resetForce;
            mesh2d.f = obj.scriptedForceInjection(mesh2d,animationScripts,frame,h);

            %gravity in y direction
            mesh2d.f(2:2:end) = mesh2d.f(2:2:end) + mesh2d.mass(2:2:end)*obj.Gravity;
            mesh2d.f = mesh2d.f  + obj.CustomForce; %+ cache.elasticForces

            %Build layers
            [rigidBodySetsPerLayer, numRigids] = buildLayers(obj,mesh2d, cache, h);
            isRigidElement = rigidBodySetsPerLayer ~= 0;
            cache.layers = isRigidElement;
            
            %finding boundaries through intersection
            [isBoundaryVertex,isElasticElement,isRigidVertex,vertexRigidBody] = xpbdLayer2D.findBoundaryVertices(obj.layers, obj.sortedLayers, mesh2d, isRigidElement, rigidBodySetsPerLayer);

            % xiF = mesh2d.getPositionFormatted;
            % pointsI = xiF(isBoundaryVertex(:,3),:);
            % plot(pointsI(:,1),pointsI(:,2),'.');

            xi = mesh2d.p;

            % doing an explicit step---------------
            mesh2d.v(mesh2d.unpinnedDOFs) = mesh2d.v(mesh2d.unpinnedDOFs) + h* mesh2d.f(mesh2d.unpinnedDOFs)./mesh2d.mass(mesh2d.unpinnedDOFs);
            xi(mesh2d.unpinnedDOFs) = xi(mesh2d.unpinnedDOFs) + h * mesh2d.v(mesh2d.unpinnedDOFs);

            %project-------------------------
            cache.xarray = xi;

            lambdai = zeros(size(mesh2d.t,1),3);


            ticNonMonitoredTime = tic;
            initialResisual = 0;
            cache.initialResidual = initialResisual;
            sumNonMonitoredTime = toc(ticNonMonitoredTime) + sumNonMonitoredTime;
            %pinned vertices are hard constraints so we ignore them in the
            %residual
            residualArray = [cache.initialResidual]; 
            contactNormals = cat(1,cInfo.normal);
            pointColliders = cat(1,cInfo.pointCollider);
            contactVertexID = cat(1,cInfo.vertexID);

            %solve  xpbdLayer2D.layerSolve or wrapLayerSolve_mex
            [xi, sumNonMonitoredTime, residualArray, obj.runUntilResSmallerThan] = xpbdLayer2D.layerSolve(isElasticElement,...
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
                                                                                                obj.Gravity, ...
                                                                                                mesh2d.perElementXPBDbeta, ...
                                                                                                obj.runOnGPU,...
                                                                                                obj.StopAtImprovementPercent, ...
                                                                                                mesh2d.elA, ...
                                                                                                obj.useRigidConstacts);
            cache.residualArray = xpbdLayer2D.prepResidual(settings,residualArray);
            cache.initialResidual = cache.residualArray(1);

            %update velocities and positions-------------------------
            mesh2d.v = BDF.BDF1(xi,mesh2d.prevp,h);
            % mesh2d.v = BDF.BDF2(xi,mesh2d.prevp,mesh2d.prevPrevp,h);
            % mesh2d.v = BDF.BDF3(xi,mesh2d.prevp,mesh2d.prevPrevp, mesh2d.prevPrevPrevp,h);
            mesh2d.p = xi;

            mesh2d.v = xpbdLayer2D.contactVelocityUpdate(contactNormals,contactVertexID, mesh2d.v, h, obj.Gravity);

            td.integrateForces = toc( ticIntegrateForces ) - sumNonMonitoredTime;
            td.residual = sum(residualArray(:,end));
            td.residualPreSolve = sum(residualArray(:,1));
        end
    end
end