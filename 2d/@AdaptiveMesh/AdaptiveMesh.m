classdef AdaptiveMesh < Mesh
    properties
        AdaptiveM           % mass matrix with rigid body support
        
        RigidBodies         % array of rigid bodies
        
        ElasticDOFs         % list of elastic DOFs in the 2 * N state vector
        ElasticInds         % indices of elastic vertices

        ElasticTriInds      % indices of elastic elements (duplication ??)
        
        RigidificationValues               % per triangle history of E dot values or other metrics, leftmost being most recent
        
        EnableAutomaticRigidification      % setting to enable rigidification
        EnableAutomaticElastification      % setting to enable elastification
        DebugEdotNorms
        VertexDisp      % local coordinates of rigid dofs in the rigid bodies
        isTriElastic    % list of #triangles booleans representing if the element is elastic or rigidified.
        ActiveElasticDOFs   % list of elastic degrees of freedom
        ActiveDofsCorrespondingID %list of the activeDofs positions for easier indexing
        RigidificationDifference %use with difference plotting, shows elements that were either implicitely rigidified or removed from rigidification (causing hinges)
        
        rigidIDbyVert %output from mexRigidificatior that gives us the rigid body to which a vertex belongs
        rigidIDbyElement
    end
    methods(Static = true)
        function [R, com, angle, inertia, rigidMasses, VertexDisp, comdot, angularVelocity]=makeRigidsFromConnectivity(numRigid, p, v, rigidIDbyVert, particlemasses)
            %USED ONLY FOR THE MULTILAYER SOLVES
            %requires mesh2d.rigidIDbyVert and rigidIDbyElement to be defined
            % [com, comdot, mass, rotMass, omega, vertexDisp] = mexRigidBodyProperties2D(numRigid, rigidIDbyVert, p, v, particlemasses);
            [com, comdot, rigidMasses, inertia, angularVelocity, VertexDisp] = RigidBody.matComputeRigidBodyProperties2D(numRigid, rigidIDbyVert, p, v, particlemasses);
            com = reshape(com',1,2,[]);
            comdot = reshape(comdot',1,2,[]);
            
            angle = zeros(1,1,numRigid);
            R = repmat(eye(2,2),[1,1,numRigid]);
        end

        [activeDOFs, ActiveElasticDOFs, ActiveDofsCorrespondingID] = computeActiveDOFs(ElasticDOFs, unpinnedDOFs, totalDOFs, RigidBodies);
    end
    methods
        function obj = AdaptiveMesh(p, t, attributes, materials)
            % ADAPTIVEMESH Extension of Mesh that supports adaptive
            % rigidification. 
            
            if nargin == 1
                % This is as ugly as it gets but I don't see any other way
                % without constructor overloading
                t = 0;
                attributes = [];
                materials = [];
            end
            
            obj@Mesh(p, t, attributes, materials);
            
            if isa(p, 'AdaptiveMesh')
                % copy constructor
                fns = properties(p);
                for i = 1:numel(fns)
                    obj.(fns{i}) = p.(fns{i});
                end
                return;
            end
            obj.VertexDisp = zeros(obj.N,2);
            obj.RigidBodies = RigidBody.empty;
            obj.ElasticDOFs = 1:obj.N * 2;
            obj.ElasticInds = 1:obj.N;
            obj.ActiveElasticDOFs = 1:obj.N*2;
            obj.rigidIDbyVert = zeros(obj.N,1);
            
            obj.ElasticTriInds = 1:size(t, 1);
            obj.AdaptiveM = obj.M;
            
            obj.EnableAutomaticRigidification = 1;
            obj.EnableAutomaticElastification = 1;
            obj.isTriElastic = true(size(obj.t,1),1);
            obj.updateRigidState(); % this should not be necessary unless something is not initialized properly but I can't find what is not!
        end
        
        
        function makeRigids(mesh2d,p,v,rigidElements)%USED ONLY FOR THE MULTILAYER SOLVES
            if numel(rigidElements)==0
                mesh2d.RigidBodies(1:end) = [];
                mesh2d.rigidIDbyVert = zeros(mesh2d.N,1);
                mesh2d.isTriElastic = true(size(mesh2d.isTriElastic));
                mesh2d.rigidIDbyVert = zeros(size(mesh2d.rigidIDbyVert));
                mesh2d.ElasticTriInds = find(mesh2d.isTriElastic)';
                mesh2d.ElasticInds = 1:mesh2d.N;
                inds2 = mesh2d.ElasticInds * 2;
                mesh2d.ElasticDOFs = reshape([inds2 - 1; inds2], 1, []);
                [mesh2d.activeDOFs, mesh2d.ActiveElasticDOFs, mesh2d.ActiveDofsCorrespondingID] = computeActiveDOFs(mesh2d.ElasticDOFs, mesh2d.unpinnedDOFs, mesh2d.N*2, mesh2d.RigidBodies);
                return;
            end
            if ~isempty(mesh2d.animationInds)                
                pinned = mesh2d.pinned;
                pinned(mesh2d.animationInds) = true;
                pinned = logical(pinned);
                pinnedVertPerTri = pinned(mesh2d.t);
                stablePinnedTri = find(sum(pinnedVertPerTri,2) >= 2);
                isStableTri = false(size(mesh2d.t,1),1);
                isStableTri(stablePinnedTri) = true;
            else
                pinned = mesh2d.pinned;
                pinned = logical(pinned);
                stablePinnedTri = mesh2d.stablePinnedTri;
                isStableTri = false(size(mesh2d.t,1),1);
                isStableTri(stablePinnedTri) = true;
            end
                
            rigidElementLogical = false(size(mesh2d.t,1),1);
            rigidElementLogical(rigidElements) = true;
            toElastify= false(size(rigidElementLogical));

            %for now forces the pinned tris to be elastic
            toElastify(mesh2d.pinnedTris) = true;

            [numRigid, mesh2d.rigidIDbyElement, mesh2d.rigidIDbyVert, isVertElastic, isVertBoundary, mesh2d.isTriElastic, isComponentPinned] = mexRigidConnectedComponents2D(rigidElementLogical, toElastify, mesh2d.AdjagencyMatrix, mesh2d.t, mesh2d.N , pinned, mesh2d.valence, stablePinnedTri, isStableTri);
           
            assert(sum(isVertElastic) == sum(mesh2d.rigidIDbyVert==0));
            assert(sum(mesh2d.isTriElastic) == sum(mesh2d.rigidIDbyElement<=0));
            assert(all(mesh2d.rigidIDbyElement~=0));
            
            mesh2d.ElasticTriInds = find(mesh2d.isTriElastic)';
            mesh2d.ElasticInds = find(isVertElastic)';
            inds2 = mesh2d.ElasticInds * 2;
            mesh2d.ElasticDOFs = reshape([inds2 - 1; inds2], 1, []);
            
            rigidVerts = find(~isVertElastic);

            [com, comdot, mass, rotMass, omega, vertexDisp] = RigidBody.matComputeRigidBodyProperties2D(numRigid, mesh2d.rigidIDbyVert, p, v, particlemasses);
            mesh2d.VertexDisp = [vertexDisp(rigidVerts*2-1),vertexDisp(rigidVerts*2)];
            
            mesh2d.RigidBodies = arrayfun(@RigidBody,1:numRigid);
            for i = 1:numRigid
                mesh2d.RigidBodies(i).TriInds = find(mesh2d.rigidIDbyElement == i)';
                mesh2d.RigidBodies(i).Indices = find(mesh2d.rigidIDbyVert == i)';
                mesh2d.RigidBodies(i).Position = com(:,i);
                mesh2d.RigidBodies(i).disgustingCom = com(:,i);
                mesh2d.RigidBodies(i).Velocity = comdot(:,i);
                mesh2d.RigidBodies(i).Mass = mass(i);
                mesh2d.RigidBodies(i).Inertia = rotMass(i);
                mesh2d.RigidBodies(i).disgustingI = mesh2d.RigidBodies(i).Inertia;
                mesh2d.RigidBodies(i).AngularVelocity = omega(i);
                mesh2d.RigidBodies(i).disgustingAngularVelocity = omega(i);
                mesh2d.RigidBodies(i).Angle = 0;
                mesh2d.RigidBodies(i).disgustingAngle = 0;
                inds2 = mesh2d.RigidBodies(i).Indices * 2;
                mesh2d.RigidBodies(i).DOFs = reshape([inds2 - 1; inds2], 1, []);
                mesh2d.RigidBodies(i).isPinned = isComponentPinned(i);
                mesh2d.RigidBodies(i).Force = [0;0];
                mesh2d.RigidBodies(i).Torque = 0;
                mesh2d.RigidBodies(i).disgustingR = eye(2);
                mesh2d.RigidBodies(i).disgustingBoundaryNum = 1;
                mesh2d.RigidBodies(i).computeri(p);
            end
            [mesh2d.activeDOFs, mesh2d.ActiveElasticDOFs, mesh2d.ActiveDofsCorrespondingID] = AdaptiveMesh.computeActiveDOFs(mesh2d.ElasticDOFs, unpinnedDOFs, numel(p), mesh2d.RigidBodies);
        end

        function clone = clone(obj)
            clone = AdaptiveMesh(obj);
        end
        
        function M = getM(obj)
            M = obj.AdaptiveM;
        end
        
        function linearMomentum = getAdaptiveLinearMomentum( mesh )
            % GETADAPTIVELINEARMOMENTUM For comparison with the easy
            % computation
            linearMomentum = sum( reshape(mesh.mass(mesh.ElasticDOFs).*mesh.v(mesh.ElasticDOFs), 2, []), 2 );
            for i = 1:numel(mesh.RigidBodies)
                linearMomentum = linearMomentum + mesh.RigidBodies(i).Mass * mesh.RigidBodies(i).Velocity;
            end
        end
        
        function labelElasticTriangles( mesh )
             pt = reshape( mesh.p, 2, mesh.N );
             for i=mesh.ElasticTriInds
                 pos = sum( pt( :, mesh.el(i).t ), 2 ) / 3;
                 text( pos(1), pos(2), ""+i );
             end        
        end


        
        prepare(obj, infMassForPinned);
        
        addForce(mesh, force);
        
        applyAcceleration(obj, acc);
        
        removeRigidBody(obj, body);
        updateParticles(mesh, h, deltav);
        
        [adaptB, gamma] = computeAdaptB(obj);
        
        v = getCurrentVelocity(mesh);
        
        f = getCurrentForce(mesh);
        
        updateRigidState(mesh);
        
        computeRigidForce(mesh,h);
    end
end

