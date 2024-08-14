classdef RigidBody < handle
    properties
        Mesh            % Mesh of this RigidBody
        
        Position        % rigid position
        Velocity        % rigid body velocity
        Angle           % rigid orientation
        AngularVelocity % rigid angular velocity

        Force           % rigid force accumulator
        Torque          % rigid torque accumulator

        Mass            % mass of rigid
        Inertia         % inertia of the rigid body
        ri              % per vertex distance vectors when the rigid body is build N by 2
        vertex2ri       % map from vertex id to ri entry

        DOFs            % list of rigid DOFs in the 2 * N state vector
        Indices         % indices of rigid vertices 
        TriInds         % rigid triangles indices in this body
        isPinned
        disgustingCom
        disgustingR
        disgustingAngle
        disgustingAngularVelocity
        disgustingBoundaryNum
        disgustingI
    end
    
    methods
        function obj = RigidBody(mesh)
            %RIGIDBODY makes a rigid body to be part of an AdaptiveMesh.
            %   Empty rigid bodies should not exist (rigid bodies with no
            %   triangles)
            if nargin == 1
                obj.Mesh = mesh;
            end
            
            obj.Position = zeros(2, 1);
            obj.Velocity = zeros(2, 1);
            obj.Angle = zeros(1, 1);
            obj.AngularVelocity = zeros(1, 1);

            obj.Force = zeros(2, 1);
            obj.Torque = zeros(1, 1);

            obj.Mass = 0;
            obj.Inertia = 0; 
            
            obj.DOFs = [];
            obj.Indices = [];
            obj.TriInds = [];  
            obj.isPinned = false;
        end

        
        function clone = clone(obj, mesh)
            clone = RigidBody(mesh);
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
            body.Velocity = body.Velocity + deltav(1:2);
            body.AngularVelocity = body.AngularVelocity + deltav(3);            
            body.Position = body.Position + h * body.Velocity;
            body.Angle = body.Angle + h * body.AngularVelocity;
            body.setDOFsFromRigid(h);
        end

        function [p,v,Velocity,AngularVelocity,Position,Angle] = peekPosition(body, h, deltav)
            Velocity = body.Velocity + deltav(1:2);
            AngularVelocity = body.AngularVelocity + deltav(3);            
            Position = body.Position + h * body.Velocity;
            Angle = body.Angle + h * body.AngularVelocity;
            [p, v] = body.peekDOFsFromRigid(Angle, Position, h);
        end
        
        function [com, R, angle, Inertia, angularv] = peekCOMRot(body, h, deltav)
            v = body.Velocity + deltav(1:2);
            angularv = body.AngularVelocity + deltav(3);
            com = body.Position + h * v;
            angle = body.Angle + h * angularv;
            Inertia = body.Inertia;
            R = angle2rotm2D(angle);
        end
        
        function pos = queryAllVertexPositions(body, rotation, centerOfMass)
            assert(size(centerOfMass,1)==1);
            pos = RigidBody.relativeVertexPosition(body.ri, rotation, centerOfMass);
        end

        function pos = queryVertexPosition(body, rotation, centerOfMass, vertex)
            assert(size(centerOfMass,1)==1);
            vertexPosr = body.vertex2ri(vertex);
            r = body.ri(vertexPosr,:);
            pos = RigidBody.relativeVertexPosition(r, rotation, centerOfMass);
        end

        function f = computeForceFromRigid(body)
            t = body.Torque;
            r1 = (body.Mesh.p(body.Indices * 2 - 1) - body.Position(1));
            r2 = (body.Mesh.p(body.Indices * 2) - body.Position(2));
            r = (r1.^2 + r2.^2).^0.5;
            force = t./r;
            f = [-force.*(r2./r) + body.Force(1),force.*(r2./r)+ body.Force(2)]./body.Mesh.mass(body.Indices);
        end
        
        function computeri(body,p)
            %must be called after Indices and Position are filled.
            body.vertex2ri = dictionary(body.Indices,1:numel(body.Indices));
            indsp = body.Indices*2;
            body.ri = [p(indsp-1),p(indsp)] - body.Position';
        end

        % public function prototypes
        computeRigidForce(body,h);
        updateBody(obj, updateMesh);
        setDOFsFromRigid(body, h);
        [p,v] = peekDOFsFromRigid(body, Angle, Position, h);
        [Force,Torque] = peekRigidForce(body,f);

    end

    methods(Static = true)

        function [R, COM, angle] = rigidSymplecticImperativeStep(COM, angle, h, COMDOT, angularv)   
            COM = COM + h .* COMDOT;
            angle = angle + h .* angularv;
            R = angle2rotm2D(angle);
        end
        
        function pos = relativeVertexPosition(ri, rotation, centerOfMass)
            assert(size(centerOfMass,1)==1);
            pos = (rotation*ri')' + centerOfMass;
        end

        function pos = pageRelativeVertexPosition(ri,rotation,centerOfMass)
            assert(size(centerOfMass,1)==1);
            pos = pagetranspose(pagemtimes(rotation,'none',ri,'transpose')) + centerOfMass;
        end

        function [com, comdot, rigidMass, inertia, omega, vertexDisp]=matComputeRigidBodyProperties2D(numRigid, rigidIDbyVert, p_in, pdot_in, massPerDOFs)
            %COULD THE PROBLEM BE SOMETHING ABOUT the mass being split
            %between vertices?
            numVert = numel(rigidIDbyVert);
            massPerVertex = massPerDOFs(1:2:end);
            p = formatPositions2D(p_in);
            pdot = formatPositions2D(pdot_in);

            isRigidVert = rigidIDbyVert~=0;
            rigidIDbyRigidVert = rigidIDbyVert(isRigidVert);
            rigidMassAccumulator = accumarray(rigidIDbyRigidVert,massPerVertex(isRigidVert),[numRigid,1]);
            rigidMass = reshape(rigidMassAccumulator,1,1,[]);
            massPos = p.*massPerVertex;
            massVel = pdot.*massPerVertex;

            sumMassPos = zeros(numRigid,2);
            sumMassPos(:,1) = accumarray(rigidIDbyRigidVert,massPos(isRigidVert,1),[numRigid,1]);
            sumMassPos(:,2) = accumarray(rigidIDbyRigidVert,massPos(isRigidVert,2),[numRigid,1]);

            sumMassVelocity = zeros(numRigid,2);
            sumMassVelocity(:,1) = accumarray(rigidIDbyRigidVert,massVel(isRigidVert,1),[numRigid,1]);
            sumMassVelocity(:,2) = accumarray(rigidIDbyRigidVert,massVel(isRigidVert,2),[numRigid,1]);
            
            com = sumMassPos./rigidMassAccumulator;
            comdot = sumMassVelocity./rigidMassAccumulator;
            
            vertexDisp = zeros(numVert,2);
            velDisp = zeros(numVert,2);

            vertexDisp(isRigidVert,:) = p(isRigidVert,:) - com(rigidIDbyRigidVert,:);
            velDisp(isRigidVert,:) = pdot(isRigidVert,:) - comdot(rigidIDbyRigidVert,:);

            %formula from http://number-none.com/blow/inertia/deriving_i.html
            rotMassComponent = massPerVertex.*dot(vertexDisp,vertexDisp,2);
            rotMassAccumulator = accumarray(rigidIDbyRigidVert,rotMassComponent(isRigidVert),[numRigid,1]);
            inertia = reshape(rotMassAccumulator,1,1,numRigid);
            
            crossPosVel = cross2D(vertexDisp,velDisp);
            rotationComponent = massPerVertex.*crossPosVel;
            omegaAccumulator = accumarray(rigidIDbyRigidVert,rotationComponent(isRigidVert),[numRigid,1]);
            omegaMomentumless = reshape(omegaAccumulator,1,1,numRigid);
            
            % preserving angular momentum
            omega=omegaMomentumless./inertia;
        end
    end
end

