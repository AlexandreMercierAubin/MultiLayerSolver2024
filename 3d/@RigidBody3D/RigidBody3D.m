classdef RigidBody3D < handle
    properties
        Mesh            % Mesh of this RigidBody
        
        Position        % rigid position
        Velocity        % rigid body velocity
        Rotation        % rigid orientation
        AngularVelocity % rigid angular velocity

        Force           % rigid force accumulator
        Torque          % rigid torque accumulator

        Mass            % mass of rigid
        Inertia         % inertia of the rigid body
        Inertia0        % body inertia at identity rotation

        DOFs            % list of rigid DOFs in the 3 * N state vector
        Indices         % indices of rigid vertices 
        TetInds         % rigid tetrahedron indices in this body
        isPinned
        isAnimated
        disgustingI %hack for xpbd
        disgustingR
        disgustingq
        disgustingCom
        disgustingDv
        r
    end
    
    methods

        function diagonalizeInertia(obj)
            [R,eigenvalueD] = eig(obj.Inertia);
            obj.Inertia = eigenvalueD;
            obj.Rotation = R;
        end

        function obj = RigidBody3D(mesh)
            %RIGIDBODY makes a rigid body to be part of an AdaptiveMesh.
            %   Empty rigid bodies should not exist (rigid bodies with no
            %   tet)
            obj.Mesh = mesh;
            obj.Position = zeros(3, 1);
            obj.Velocity = zeros(3, 1);
            obj.Rotation = eye(3);

            obj.Force = zeros(3, 1);
            obj.Torque = zeros(3, 1);

            obj.Mass = 0;
            obj.Inertia = zeros(3,3); 
            obj.Inertia0 = zeros(3,3); 
            
            obj.DOFs = [];
            obj.Indices = [];
            obj.TetInds = [];
            obj.AngularVelocity = zeros(3,1);
            obj.isPinned = false;
        end
        
        function clone = clone(obj, mesh)
            clone = RigidBody3D(mesh);
            fns = properties(obj);
            for i = 1:numel(fns)
                clone.(fns{i}) = obj.(fns{i});
            end
            clone.Mesh = mesh;
        end
        
        function addVertices(obj, indices, updateMesh)
            if nargin < 3
                updateMesh = 1;
            end
            
            obj.Indices = unique([obj.Indices, indices]);
            obj.updateBody(updateMesh);
        end
        
        function removeVertices(obj, indices)
            obj.Indices(ismember(indices, obj.Indices)) = [];
            obj.updateBody();
        end
        
        removeTriangles(obj, tris, updateMesh, settings);
        
        function merge(obj, body, updateMesh)
            if nargin < 3
                updateMesh = 1;
            end
            obj.Mesh.RigidBodies(obj.Mesh.RigidBodies == body) = [];
            obj.Indices = unique([obj.Indices, body.Indices]);
            obj.updateBody(updateMesh);
        end
        
        function updatePosition(body, h, deltav)
            % updatePosition Performs the symplectic update of velocity and
            % the position for both linear and angular velocity
            body.Velocity = body.Velocity + deltav(1:3);
            body.Position = body.Position + h * body.Velocity;
            if body.isAnimated
                %prevent rotation for scripted animations
                body.AngularVelocity = body.AngularVelocity;
            else
                body.AngularVelocity = body.AngularVelocity + deltav(4:6);
            end
            body.Rotation = expRodrigues( body.AngularVelocity, h ) * body.Rotation;
            body.Inertia = body.Rotation * body.Inertia0 * body.Rotation'; 
            body.setDOFsFromRigid(h);
        end

        function [pList, vList,inds, com, angularv,Inertia] = peekPosition(body, h, deltav, mesh3Dp)
            % This is just to see the positions, it does not update rigid
            % body properties
            [com, R, Inertia,angularv] = peekCOMRot(body, h, deltav);
            [pList, vList, inds] = body.peekDOFsFromRigid(h, R, bodyp, mesh3Dp);
        end

        function [com, R, Inertia, angularv] = peekCOMRot(body, h, deltav)
            com = body.Position + h * deltav(1:3);
            if body.isAnimated
                %prevent rotation for scripted animations
                angularv = body.AngularVelocity;
            else
                angularv = body.AngularVelocity + deltav(4:6);
            end
            R = expRodrigues( angularv, h ) * body.Rotation;
            Inertia = R * body.Inertia0 * R';
        end

       
        %TODO: move this to the rigidbody update thingy so it is only computed
        %when rigid bodies are changed. This function will only extract the
        %values of alpha0
        function alpha0 = getAlpha0(body, gamma)
            ids = body.Indices * 3;
            alpha0s = zeros(size(ids,2)*3,1);
            alpha0s(3:3:end) = body.Mesh.alpha0(ids);
            alpha0s(2:3:end) = body.Mesh.alpha0(ids-1);
            alpha0s(1:3:end) = body.Mesh.alpha0(ids-2);
            croppedGamma = zeros(size(ids,2)*3,3);
            croppedGamma(3:3:end) = gamma(ids,:);
            croppedGamma(2:3:end) = gamma(ids-1,:);
            croppedGamma(1:3:end) = gamma(ids-2,:);
            alpha0 = zeros(3,1);
            
            for j = 1:1:size(gamma,2)
                alpha0(j) = sum(alpha0s.*croppedGamma(:,j));
            end
        end

        function pos = queryAllVertexPositions(obj, initialp, initialCenterOfMass, rotation, centerOfMass) %material pos should be the position at construction of the vertex
            xList = obj.Mesh.formatPositions(initialp(obj.DOFs));
            rList = (xList - initialCenterOfMass)';
            pos = (rotation*rList)' + centerOfMass;
        end

        

        function pos = queryAllVertexPositionsFromQuat(obj, initialp, initialCenterOfMass, quat, centerOfMass) %material pos should be the position at construction of the vertex
            xList = obj.Mesh.formatPositions(initialp(obj.DOFs));
            rList = (xList - initialCenterOfMass);
%             rotatedp = rotatepoint(quat,rList);
            rotatedp = mexQuatRot(rList, quat);
            pos = rotatedp + centerOfMass;
        end

        computeRigidForce( body, h);
        [Force, Torque] = peekRigidForce(body, h, AngularVelocity, RigidBodyPosition, Inertia, meshf, meshp);
        updateBody( obj, updateMesh );  
        setDOFsFromRigid( body, h );
    end
    methods(Static)
        
        function [R, COM, Inertia] = rigidSymplecticImperativeStep(COM, R, h, COMDOT, angularv, Inertia0)   
            COM = COM + h .* COMDOT;
            deltaR = expRodrigues( angularv, h );
            R =  pagemtimes(deltaR,R);
            Inertia = pagemtimes(pagemtimes(R, Inertia0), 'none', R,'transpose');
        end
        
        function pos = queryVertexPosition(prevp, rotation, initialCenterOfMass, centerOfMass) %material pos should be the position at construction of the vertex
            r = (prevp - initialCenterOfMass)';
            pos = (rotation*r)' + centerOfMass;
        end

        function pos = queryVertexPositionFromQuat(prevp, quat, initialCenterOfMass, centerOfMass) %material pos should be the position at construction of the vertex
            r = (prevp - initialCenterOfMass);
%             rotatedp = rotatepoint(quat,r);
            rotatedp = mexQuatRot(r, quat);
            pos = rotatedp + centerOfMass;
        end

        function pos = pageRelativeVertexPosition(ri,rotation,centerOfMass)
            assert(size(centerOfMass,1)==1);
            pos = pagetranspose(pagemtimes(rotation,'none',ri,'transpose')) + centerOfMass;
        end

        function [com, comdot, rigidMass, inertia, omega, vertexDisp]=matComputeRigidBodyProperties3D(numRigid, rigidIDbyVert, p_in, pdot_in, diagMass)
            numVert = numel(rigidIDbyVert);
            massPerVert = diagMass(1:3:end);
            pageMassPerVert = reshape(massPerVert,1,1,[]);
            p = formatPositions3D(p_in);
            pPage = pagePositions3D(p_in);
            pdot = formatPositions3D(pdot_in);
            pdotPage = pagePositions3D(pdot_in);

            isRigidVert = rigidIDbyVert~=0;
            rigidIDbyRigidVert = rigidIDbyVert(isRigidVert);
            rigidMassAcc = accumarray(rigidIDbyRigidVert,massPerVert(isRigidVert),[numRigid,1]);
            rigidMass = reshape(rigidMassAcc,1,1,[]);

            massPos = p.*massPerVert;
            massVel = pdot.*massPerVert;

            com=zeros(1,3,numRigid);
            com(1,1,:) = accumarray(rigidIDbyRigidVert,massPos(isRigidVert,1),[numRigid,1]);
            com(1,2,:) = accumarray(rigidIDbyRigidVert,massPos(isRigidVert,2),[numRigid,1]);
            com(1,3,:) = accumarray(rigidIDbyRigidVert,massPos(isRigidVert,3),[numRigid,1]);

            comdot=zeros(1,3,numRigid);
            comdot(1,1,:) = accumarray(rigidIDbyRigidVert,massVel(isRigidVert,1),[numRigid,1]);
            comdot(1,2,:) = accumarray(rigidIDbyRigidVert,massVel(isRigidVert,2),[numRigid,1]);
            comdot(1,3,:) = accumarray(rigidIDbyRigidVert,massVel(isRigidVert,3),[numRigid,1]);
            
            com = com./rigidMass;
            comdot = comdot./rigidMass;
            
            % easier to pick up display positions, J and omega in a second pass
            % compute omega so as to conserve angular momentum
            vertexDisp = zeros(1,3,numVert);
            velDisp = zeros(1,3,numVert);

            vertexDisp(:,:,isRigidVert) = pPage(:,:,isRigidVert) - com(:,:,rigidIDbyRigidVert);
            velDisp(:,:,isRigidVert) = pdotPage(:,:,isRigidVert) - comdot(:,:,rigidIDbyRigidVert);

            inertia = RigidBody3D.buildInertiaTensor(vertexDisp, numRigid, pageMassPerVert, rigidIDbyRigidVert, isRigidVert);

            componentAngularMomentum = pageMassPerVert.*cross(vertexDisp,velDisp,2);
            angularMomentum = zeros(3,1,numRigid);
            angularMomentum(1,1,:) = accumarray(rigidIDbyRigidVert, reshape(componentAngularMomentum(1,1,isRigidVert),[],1,1));
            angularMomentum(2,1,:)  = accumarray(rigidIDbyRigidVert, reshape(componentAngularMomentum(1,2,isRigidVert),[],1,1));
            angularMomentum(3,1,:)  = accumarray(rigidIDbyRigidVert, reshape(componentAngularMomentum(1,3,isRigidVert),[],1,1));
            
            % a final pass to fix up our mass weighted computation of omega (i.e., preserving angular momentum)
            angularVels = pagemldivide(inertia,angularMomentum);
            omega=pagetranspose(angularVels);
        end

        function inertia = buildInertiaTensor(vertexDisp, numRigid, perVertMassArray, rigidIDbyRigidVert, isRigidVert)
            rhat = crossProductMatrix(vertexDisp);
            rotMassComponent = -perVertMassArray .* pagemtimes(rhat,'none',rhat,'none');

            inertia = zeros(3,3,numRigid);%TODO: there must be a better way to use that function here without a for loop as the indexing is the same
            for inertiaIndexx = 1:3
                for inertiaIndexy = 1:3
                    inertia(inertiaIndexx,inertiaIndexy,:) = accumarray(rigidIDbyRigidVert,reshape(rotMassComponent(inertiaIndexx,inertiaIndexy,isRigidVert),[],1,1));
                end
            end

            % And also the omega cross J omega term...
            % L = Inertia*AngularVelocity;
            % Lhat = crossProductMatrix(L)*h;
            % Jmod = Inertia - h*h*Lhat*diag(1./[body.Mass,body.Mass,body.Mass])*Lhat;
            % Torque = Torque - cross( AngularVelocity, Jmod * AngularVelocity );
        end
    end
end

