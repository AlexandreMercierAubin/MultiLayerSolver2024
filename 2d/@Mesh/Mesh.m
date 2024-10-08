classdef Mesh < handle
    properties
        N   % number of DOFs
        p   % node positions
        p0  % initial node positions
        prevp
        prevPrevp
        prevPrevPrevp
        v   % node velocities
        v0  % initial node velocities
        f   % force accumulator for nodes
        
        t      % triangles (for making elements)
        tint32 % indices as integers
        elementDOFs %just a different representation of t

        el     % elements
        boundaryEdges  % cell array of boundaries (for collision detection)
        % Every mesh has at least one set of boundary edges
        % Each edge is a line, first column is the index of the first point
        % second column is the index of the second point
        % NOTE: order is important!!!  Treating edges as directed, the mesh
        % is to the left andt the outside is to the right, e.g., CCW order
        % for convex polygons.
        boundaryEdgeSets % which sets of boundaries can be checked with each other for collision
        
        Graph                   % Triangle adjacency graph (across edges)
        AdjagencyMatrix         % Adjacency matrix used to build the graph (not used otherwise?)
        TrianglesPerParticle   % cell array of the elements adjacent to a particle
        
        pinned       % flags for pinned indices
        pinnedInds   % indices of pinned vertices
        pinnedDOFs   % pinned indices in the 2 * N state vector
        unpinnedDOFs % unpinned indices in the 2 * N state vector
        objectDOFs %cell array of lists that containts the dofs of each individual objects
        objectDOFsPinned %same, but with pinned DOFs removed
        pinnedTris   % pinned triangles (triangles with at least a pinned particle)
        isTriPinned  % logical array for pinned triangles
        stablePinnedTri %tri with at least 2 vertex pinned
        isStableTri % same, but in logical array format
        isVertexPinned  % logical array for vertices of triangles
        cacheMD5 % MD5 solely use to detect change in merged meshes scene creation
        cacheMergeMaterials %solely use to detect change in merged meshes scene creation
        animationDOFs
        animationInds

        perElementXPBDalphaMatrix
        perElementXPBDalphaInvMatrix
        perElementXPBDbeta
        
        activeDOFs
        ActiveBRows % indices of rows in B which are elastic (active), i.e., 1 2 3 4 if only the first element
        
        mass    % diag of mass matrix (likely more convenient?)
        M       % sparse mass matrix, 2 * N by 2 * N    
        Mii
        enforcedMassIndices %logical array of indices to set mass instead of lumping. Needs to be set before the material update
        enforcedMass % the mass of the particles forced to be a given value
        
        alpha0mass     % diag of lumped Rayleigh alpha0 multiplied with lumped mass
        Md             % sparse matrix version of alpha0mass   
        Mdii
        
        B              % kinematic relationship between vertices and element deformation gradients F
        DmInv          % inverse element frame
        Bii

        %TODO: THIS NEVER WORKED, and should probably be stripped out of
        %the code along with all the other mess associated with it (i.e.,
        %code in plotMesh, etc.
        lagrangeMults  % lambdas for the compliant feedback constraint version of elasticity

        %Cf   % symbolic d2psid2F function - deprectated use mCSTVK
        
        elBm        % quick access to the shape matrix for all elements (2x2xnumel)
        elA         % quick per element area access for C computation
        elMu        % quick per element lame access for C computation
        elGamma     % for stable neohookean constraints's compliance
        elLambda    % quick per element lambda access for C computation
        elAlpha1    % quick per element alpha 1 acces
        
        bigAlpha1   % sparse matrix for combining with full C to build Kd
        
        RenderPatch        % plot of the elements
        BoundaryLines      % plot of the boundaries

        % Note: this resuablefigures member may not the cleanest way to
        % cache figure handles. Pick a name for a figure, can make it a 
        % member of the reusableFigures struct.  Use a combination of 
        % isfield(reusableFigs,'name') isvalid(reusableFigs.name) and
        % delete(reusableFigs.name) to manage your figures.   
        reusableFigs = struct('none',0);

        materials           % material per triangle
        materialIndex       % triangle index to material
        
        valence %valence of the vertices
        edges
        edgeRest
        plotOffset = [0,0];

        elementColor %graph coloring
        elementRGB
        numColors;
    end
    
    methods
        function obj = Mesh(p, t, attributes, materials, enforcedMassIndices, enforcedMass)
            % MESH Constructs a mesh from given points and triangles.
            %   Prepares for simulation with provided materials.
            %
            %   Note that if you use meshes saved in a .mat file, you will
            %   need to regenerate the mesh if new fields get added
            %   
            if nargin == 2 || (nargin >= 3 && size(attributes,2) <= 0)
                attributes = ones(size(t,1),1);
            end
            
            if nargin <= 3 || (nargin >= 4 && size(materials,2) <= 0)
                defaultMaterial = TriangleMaterial();
            else
                defaultMaterial = materials;
            end
            
            if nargin >= 2 && size(attributes,2) > size(materials,2)
                ids = setdiff(1:size(attributes,2),1:size(materials,2));
                defaultMaterial(ids) = TriangleMaterial();
            end

            if nargin < 5
                enforcedMassIndices = false(size(p(:)));
                enforcedMass = -1;
            end
            
            if isa(p, 'Mesh')
               % copy constructor
                fns = properties(p);
                for i = 1:numel(fns)
                    obj.(fns{i}) = p.(fns{i});
                end
                return;
            end
            N = size(p, 1);
            obj.N = N;                 % number of nodes
            obj.objectDOFs = {1:N*2};
            obj.objectDOFsPinned = obj.objectDOFs;
            obj.p = reshape(p', N * 2, 1); % position state
            obj.prevp = obj.p;
            obj.prevPrevp = obj.p;
            obj.prevPrevPrevp = obj.p;
            obj.p0 = obj.p;           % initial positions
            obj.v = zeros(N * 2, 1);      % velocity state
            obj.v0 = zeros(N * 2, 1);     % initial velocities
            obj.f = zeros(N * 2, 1);      % force accumulator
            obj.enforcedMassIndices = enforcedMassIndices;
            obj.enforcedMass = enforcedMass;

            obj.t = int32(t);                 % triangles
            obj.elementDOFs = [obj.t(:,1)*2-1,obj.t(:,1)*2,obj.t(:,2)*2-1,obj.t(:,2)*2,obj.t(:,3)*2-1,obj.t(:,3)*2];
            
            obj.el = obj.makeElements(p);
            obj.boundaryEdges = obj.makeBoundaries(p);
            obj.boundaryEdgeSets = { 1:numel(obj.boundaryEdges) };
            
            obj.elementColor = naiveGraphColoring(obj.t);
            obj.numColors = max(obj.elementColor);
            obj.elementRGB = rand(obj.numColors,3);

            %computing valence
            obj.valence = zeros(obj.N,1);
            for i = 1:size(obj.t,1)
                for j = 1:size(obj.t,2)
                    vertex = obj.t(i,j);
                    obj.valence(vertex) = obj.valence(vertex) + 1;
                end
            end
            
            obj.updateMaterials(attributes,[obj.materials, defaultMaterial]);
            [obj.AdjagencyMatrix,obj.TrianglesPerParticle,obj.Graph] = elementAdjacencyMatrix(obj.t,N);
            
            obj.pinned = zeros(N, 1);   % flags pinned indices
            obj.pinnedInds = [];        % list of pinned node IDs
            obj.pinnedDOFs = [];
            obj.unpinnedDOFs = 1:N*2;
            obj.pinnedTris = [];
            obj.isTriPinned = zeros(size(obj.t,1),1);
            obj.activeDOFs = 1:2 * N;
            obj.ActiveBRows = 1:size(t, 1) * 4;

            [obj.B,obj.DmInv] = obj.computeB2D();
            obj.Bii = obj.B(:,obj.unpinnedDOFs);
            obj.lagrangeMults = zeros(size(obj.B, 1), 1);

            %updateCf(obj);
            
            obj.elBm = [ obj.el(:).Bm ];
            obj.elA = [ obj.el(:).area ];
            

            obj.RenderPatch = 0;
            obj.BoundaryLines = {};

            duplicateEdge = containers.Map('KeyType','char','ValueType','double');
            obj.edges = [];
            for i = 1: size(obj.t,1)
                vert1 = obj.t(i,1);
                vert2 = obj.t(i,2);
                vert3 = obj.t(i,3);
                if ~duplicateEdge.isKey(char([vert2,vert3])) && ~duplicateEdge.isKey(char([vert3,vert2])) 
                    obj.edges = [obj.edges; vert2, vert3];
                    duplicateEdge(char([vert3,vert2])) = true;
                end
                
                if ~duplicateEdge.isKey(char([vert1,vert2])) && ~duplicateEdge.isKey(char([vert2,vert1])) 
                    obj.edges = [obj.edges; vert1, vert2];
                    duplicateEdge(char([vert1,vert2])) = true;
                end
                
                if ~duplicateEdge.isKey(char([vert1,vert3])) && ~duplicateEdge.isKey(char([vert3,vert1])) 
                    obj.edges = [obj.edges; vert1, vert3];
                    duplicateEdge(char([vert1,vert3])) = true;
                end
            end
            edgep1 = p(obj.edges(:,1),:);
            edgep2 = p(obj.edges(:,2),:);
            obj.edgeRest = vecnorm(edgep2-edgep1,2,2);
        end
        
        function updateMaterials( obj, newAttributes, newMaterials )        
            obj.materials = newMaterials;
            obj.materialIndex = newAttributes;
            
            updateMass(obj);
            obj.M = sparse( 1:2*obj.N, 1:2*obj.N, obj.mass ); % sparse matrix form for convenience
            obj.Mii = obj.M;
            updateAlpha0(obj);
            obj.Md = sparse( 1:2*obj.N, 1:2*obj.N, obj.alpha0mass ); % sparse matrix form for convenience
            obj.Mdii = obj.Md;

            obj.elA = [ obj.el(:).area ];
            obj.elLambda = [ obj.materials(obj.materialIndex(:)).lambda ];
            obj.elAlpha1 = [ obj.materials(obj.materialIndex(:)).alpha1 ];
            obj.elMu = [ obj.materials(obj.materialIndex(:)).mu ];
            obj.elGamma = [ 1 + obj.elMu ./ obj.elLambda ];
            obj.elBm = [ obj.el(:).Bm ];

            obj.perElementXPBDalphaMatrix = zeros(3,3,size(obj.t,1));
            obj.perElementXPBDalphaInvMatrix = zeros(3,3,size(obj.t,1));
            for i = 1:size(obj.t,1)
                lam=obj.elLambda(i);
                diagVal = lam + 2*obj.elMu(i);
                matK = [diagVal,lam,0;
                        lam, diagVal,0;
                        0,0,2*obj.elMu(i)];
                obj.perElementXPBDalphaInvMatrix(:,:,i)=matK; %confusing, but the one inverted in the uninverted one in the paper
                obj.perElementXPBDbeta(:,:,i) = matK.*obj.elAlpha1(i);
            end
            obj.perElementXPBDalphaMatrix = pageinv(obj.perElementXPBDalphaInvMatrix);

            alpha1s = reshape( repmat( obj.elAlpha1, 4, 1 ), [], 1 );
            obj.bigAlpha1 = sparse( 1:numel(alpha1s), 1:numel(alpha1s), alpha1s );
        end
        
        function mergeMesh(obj, mesh) 
            % MergeMesh adds a given mesh to this mesh
            % this is not necessarily efficient in that it will recompute
            % things that were computed before, but it all precomputation
            % so perhaps we don't need to worry so much.
            
            N1 = obj.N;
            N2 = mesh.N;
            numTris1 = size(obj.t,1);
            
            obj.N = N1 + N2;
            obj.objectDOFs = {obj.objectDOFs{:}, mesh.objectDOFs{:} + (N1*2)};
            obj.objectDOFsPinned = obj.objectDOFs;
            obj.p = [ obj.p; mesh.p ];
            obj.prevp = obj.p;
            obj.prevPrevp = obj.p;
            obj.prevPrevPrevp = obj.p;
            obj.p0 = [ obj.p0; mesh.p0 ];
            obj.v = [ obj.v; mesh.v ];
            obj.v0 = [ obj.v0; mesh.v0 ];
            obj.f = [ obj.f; mesh.f ];
            
            % triangle indices of the added mesh are offset by the number
            % of particles in this mesh
            obj.t = [ obj.t; mesh.t + N1 ];  % triangles
            obj.elementDOFs = [obj.elementDOFs; mesh.elementDOFs + (N1*2)];
            obj.valence = [obj.valence; mesh.valence];
            obj.elementColor = [obj.elementColor;mesh.elementColor];
            obj.numColors = max(obj.elementColor);
            obj.elementRGB = rand(obj.numColors,3);
            % TODO: don't really need to remake them, see fix below:
            % easiest because the elements stuipdly contain t... should 
            % refactor so that the mesh simply has Bm and area and 
            % make elements is simply a private method or part of the
            % constructor.
            obj.el = obj.makeElements( reshape( obj.p0', 2, [] )');
   
            numBoundaryEdges1 = numel(obj.boundaryEdges);
            for j = 1:numel(mesh.boundaryEdges)
                obj.boundaryEdges{numBoundaryEdges1+j} = mesh.boundaryEdges{j} + N1; 
            end
            
            numBoundaryEdgeSets1 = numel(obj.boundaryEdgeSets);
            for j = 1:numel(mesh.boundaryEdgeSets)
                obj.boundaryEdgeSets{numBoundaryEdgeSets1+j} = mesh.boundaryEdgeSets{j} + numBoundaryEdges1;
            end
            
            numMat1 = numel(obj.materials);
            obj.materials = [ obj.materials, mesh.materials ];
            obj.materialIndex = [ obj.materialIndex; mesh.materialIndex + numMat1 ];
            
            obj.updateMaterials(obj.materialIndex, obj.materials);
            
            % concat the triangle adjacency matrix and graph
            obj.AdjagencyMatrix = blkdiag(obj.AdjagencyMatrix,mesh.AdjagencyMatrix);
            obj.TrianglesPerParticle = [obj.TrianglesPerParticle;mesh.TrianglesPerParticle];
                  
            obj.activeDOFs = 1:2 * obj.N;
            obj.ActiveBRows = 1:size(obj.t, 1) * 4;
            
            obj.B = blkdiag(obj.B,mesh.B);
            obj.DmInv = cat(3,obj.DmInv,mesh.DmInv);
            obj.Bii = obj.B(:,obj.unpinnedDOFs);
            obj.lagrangeMults = zeros(size(obj.B, 1), 1);

            %updateCf(obj);
            
            obj.RenderPatch = 0;
            obj.BoundaryLines = {};
            obj.pin([mesh.pinnedTris + N1]');
            obj.edges = [obj.edges;mesh.edges];
            obj.edgeRest = [obj.edgeRest;mesh.edgeRest];
        end
        
        function x = getPositionFormatted(obj)
            %Get x in a n x 2 format compatible with most gptoolbox
            %functions
            x = reshape(obj.p',2,obj.N)';
        end

        function clone = clone(obj)
            clone = Mesh(obj);
        end
        
        function [B, gamma] = getB(obj, ~)
           B = obj.B;
           gamma = speye(obj.N * 2);
        end
        
        function M = getM(mesh)
            M = mesh.M;
        end
        
        function f = getCurrentForce(mesh)
            f = mesh.f;
        end
        
        function v = getCurrentVelocity(mesh)
            v = mesh.v;
        end
        
        function linearMomentum = getLinearMomentum( mesh )
            linearMomentum = sum( reshape(mesh.mass.*mesh.v, 2, mesh.N), 2 );
        end
        
        function labelTriangles( mesh )
             pt = reshape( mesh.p, 2, mesh.N );
             for i=1:numel(mesh.el)
                 pos = sum( pt( :, mesh.el(i).t ), 2 ) / 3;
                 text( pos(1), pos(2), ""+i );
             end        
        end

        function labelVertices( mesh )
             pt = reshape( mesh.p, 2, mesh.N );
             for i=1:mesh.N
                 pos = pt(:,i);
                 text( pos(1), pos(2), ""+i );
             end        
        end

        function updateMass(mesh)
            % UPDATEMASS is used to contsruct the mass damping and is
            % typically only called on construction of the mesh
            mesh.mass = computeMass(mesh);
        end
        
        function updateAlpha0( mesh )
            % UPDATEALPHA0 is used to contsruct the mass damping and is
            % typically only called on construction of the mesh
            
            alpha0triangle = [mesh.materials(mesh.materialIndex).alpha0];
            %values are declared per triangle, but we want them per node
            alpha0vector = zeros( mesh.N*2, 1 );
            for i = 1:size(mesh.t,1)
                for node = 1:size(mesh.t,2)
                    id = mesh.t(i,node);
                    value = alpha0vector(id*2) + 1/3 * alpha0triangle(i);
                    alpha0vector(id*2) = value;
                    alpha0vector(id*2-1) = value;
                end
            end
            mesh.alpha0mass = alpha0vector .* mesh.mass ;
        end
        
        function V = getEnergy( mesh )
            F = mesh.B * mesh.p;
            F = reshape( F, 2, 2, [] );
            n = size(F, 3);
            I = [ones(1, 1, n), zeros(1, 1, n); zeros(1, 1, n), ones(1, 1, n)];
            E = 1 / 2 * (mult2x2TransposeIn3D(F, F) - I);
            E = reshape( E, 4, [] );
            psi = mesh.elMu .* sum( E.*E, 1 ) + mesh.elLambda/2 .* (E(1,:).*E(4,:)).^2;
            V = sum( mesh.elA .* psi );
        end
        
%         function updateCf(mesh)
%         % Deprectated... use mCSTVK instead
%         % Might be interesting to keep this in comments for ease of
%         % creating a new energy??
%             Fs = sym('Fs', [2, 2],'real');  
%             E = 1 / 2 * (Fs' * Fs - eye(2));
%             mus = sym('mus','real');
%             lambdas = sym('lambdas','real');
%             P = Fs * (2 * mus * E + lambdas * trace(E) * eye(2));
%             Pr = reshape(P, 4, 1);
%             Fr = reshape(Fs, 1, 4);
% 
%             Cs = sym(zeros(4, 4));
%             for j = 1:4
%                 Cs(:, j) = diff(Pr, Fr(j));
%             end
%             mesh.Cf = matlabFunction(Cs);%, mus, lambdas
%         end
        
function [B,Dminv] = computeB2D(mesh2d)
            %COMPUTEB Slow but only precomputation. Could probably be vectorized.
            %Not a huge bottleneck for gigantic meshes (adjagency matrix is worse)
            n = size(mesh2d.p, 1);
            m = size(mesh2d.t, 1);
            B = sparse(size(mesh2d.t, 1), n);
            
            Dm = zeros(2, 2, m);
            
            %p0 - p2
            Dm(1, 1, :) = mesh2d.p0(mesh2d.t(:, 1) * 2 - 1) - mesh2d.p0(mesh2d.t(:, 3) * 2 - 1);
            Dm(2, 1, :) = mesh2d.p0(mesh2d.t(:, 1) * 2) - mesh2d.p0(mesh2d.t(:, 3) * 2);
            %p1 - p2
            Dm(1, 2, :) = mesh2d.p0(mesh2d.t(:, 2) * 2 - 1) - mesh2d.p0(mesh2d.t(:, 3) * 2 - 1);
            Dm(2, 2, :) = mesh2d.p0(mesh2d.t(:, 2) * 2) - mesh2d.p0(mesh2d.t(:, 3) * 2);
            
            cells = cellfun(@inv, num2cell(Dm, [1 2]), 'UniformOutput', false);
            Dminv = cat(3, cells{:});
            
            % Vectorize this if you can...
            for triangleID = 1:size(mesh2d.t, 1)
                
                Ds = zeros(2, 2, n);
                Ds(1, 1, :) = (mesh2d.t(triangleID, 1) * 2 - 1 == 1:n) - (mesh2d.t(triangleID, 3) * 2 - 1 == 1:n);
                Ds(2, 1, :) = (mesh2d.t(triangleID, 1) * 2 == 1:n) - (mesh2d.t(triangleID, 3) * 2 == 1:n); 
                Ds(1, 2, :) = (mesh2d.t(triangleID, 2) * 2 - 1 == 1:n) - (mesh2d.t(triangleID, 3) * 2 - 1 == 1:n);
                Ds(2, 2, :) = (mesh2d.t(triangleID, 2) * 2 == 1:n) - (mesh2d.t(triangleID, 3) * 2 == 1:n);
                
                B(triangleID * 4 - 3:triangleID * 4, 1:n) = reshape(mult2x2Stack(Ds, Dminv(:, :, triangleID)), 4, n);
            end
        end

        function el = makeElements(mesh2d,p)
            % MAKEELEMENTS Creates elements from triangles.
            %   makeElements( p, t ) returns an array of structures given 2D
            %   points p (stored as rows), and triangles t defined by indices
            %   (1 indexed and stored as rows).
            %   p is #points by 2
            %   t is #elements by 3
            t = mesh2d.t;
            el = repmat(struct('t', 1), size(t, 1), 1);
            for index = 1:size(t, 1)
                i = t(index, 1);
                j = t(index, 2);
                k = t(index, 3);
        
                Dm = [p(j,:)' - p(i,:)', p(k,:)' - p(i,:)'];
        
                el(index).t = t(index,:);
                el(index).Bm = inv(Dm);
                el(index).area = 1 / 2 * abs(det(Dm));
                el(index).restLength = Dm;
            end
        end

        function boundaries = makeBoundaries(mesh2d, p)
            % makeBoundaries Creates a cell array of boundareis from triangles.
            %   makeBoundaries( p, t ) returns a cell array of boundaries, each
            %   defined by edges in connected order, computed from given 2D
            %   points p (stored as rows), and triangles t defined by indices
            %   (1 indexed and stored as rows).
            % 
            %   Assumes there is less than 100000 points in the mesh (gross yes).
            t = mesh2d.t;
            N = 100000;
            if size(p, 1) >= N
                error('triangles must have less than 100000 points');
            end
            s = java.util.HashSet();
            for i = 1:size(t, 1)
                for j = 1:3
                    a = t(i, j);
                    b = t(i, mod(j, 3) + 1);
                    hash1 = a * N + b;
                    hash2 = b * N + a;
                    if s.contains(hash2)
                        s.remove(hash2);
                    else
                        s.add(hash1);
                    end
                end
            end
            L = s.toArray;
            edges = zeros(size(L, 1), 2);
            for i = 1:size(L, 1)
                edges(i,:) = [floor(L(i) / N), mod(L(i), N)];
            end
            
            %sorting edges
            boundaries = cell(1,0);
            newEdges = [edges(1,:)];
            edges(1,:) = [];
            while ~isempty(edges)        
                for j = 1:1:numel(edges)
                    if newEdges(end,2) == edges(j,1)
                        newEdges = [newEdges; edges(j,:)];
                        edges(j,:) = [];
                        % check if we have come full circle... 
                        % if yes, save the boundary!
                        % check if there are more, and if yes, set it upt
                        if ( newEdges(end,2) == newEdges(1,1) )
                            boundaries = horzcat(boundaries, newEdges);
                            if ~isempty(edges)
                                newEdges = [edges(1,:)];
                                edges(1,:) = [];
                            end
                        end
                        break;
                    end
                end
            end
        end

        function x = formatPositions(obj,p)
            x=formatPositions2D(p);
        end

        %% public function prototypes
        prepare( obj, infMassForPinned )
        pin( obj, pinnedInds )
        updateParticles( mesh2d, h, deltav )
        applyAcceleration( obj, acc )
        resetForce( mesh2d )
        setRigidTransform( mesh2d, degrees, translation, forceCenterOfMassRotation)
        setRigidMotion( mesh2d, radPerSec, velocity, h, debugPlot)
        computeActiveDOFsInitial(mesh2d)
    end
end

