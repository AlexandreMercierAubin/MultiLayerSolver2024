classdef MeshSCD3D < ContactFinder
    % MeshSCD3D Class for mesh self collision detection.
    %   The SCD test is not general, and instead works only for the case
    %   where the mesh consists of multiple connected components, and tests
    %   the boundaries betweeen those components.
    
    properties
        quality = false; % speed instead of quality by default
        skipset = [];
        epsilon = 0.0001;
    end
    
    methods
        function obj = MeshSCD3D( frictionCoefficient, quality , skipset)
            if nargin >= 1
                obj.FrictionCoefficient = frictionCoefficient;
            end
            if nargin >= 2
                obj.quality = quality;
            end
            if nargin >=3
                obj.skipset = skipset;
            end
        end
        
        function [Jc, phi, cInfo] = findContacts( obj, mesh3D, ~, p)
            % Input mesh with N DOFS
            % Output Jc will have N columns with contact * 2  rows
            % phi will be a column vector of the size contacts
            % contactInfo is a cell array with size equal to number of contacts
            
            phi = zeros( 0, 1 );
            Jc = zeros( 0, mesh3D.N*3 );
            cInfo = contactInfo3D.empty;
            % find collisions
            collisionPairs = [];
            AABBpairs = [];
            
            for facesCellIndexObject1 = 1:numel(mesh3D.facesSets)
                if any(facesCellIndexObject1 == obj.skipset)
                    continue
                end
                facesOtherObjects = facesCellIndexObject1+1:numel(mesh3D.facesSets);

                for facesCellIndexObject2 = facesOtherObjects
                    if facesCellIndexObject1 == facesCellIndexObject2 || any(facesCellIndexObject2 == obj.skipset)
                        continue;
                    end

                    %Broad phase AABB boxes check
                    bes1Vert = 1:mesh3D.mergeMeshVerticeCount(facesCellIndexObject1);
                    bes2Vert = 1:mesh3D.mergeMeshVerticeCount(facesCellIndexObject2);
                    
                    adjustedBes1 = bes1Vert+ mesh3D.facesSetsOffset(facesCellIndexObject1);
                    adjustedBes2 = bes2Vert+ mesh3D.facesSetsOffset(facesCellIndexObject2);
                    
                    s1BoundaryPoints     = [ p(adjustedBes1*3-2), p(adjustedBes1*3-1), p(adjustedBes1*3) ];
                    s2BoundaryPoints     = [ p(adjustedBes2*3-2), p(adjustedBes2*3-1), p(adjustedBes2*3) ];
        
                    AABBpair = [min( s1BoundaryPoints,[],1), max(s1BoundaryPoints,[],1),min( s2BoundaryPoints,[],1),max(s2BoundaryPoints,[],1)];
                    [isColliding,~, ~, ~] = AABBboxesCollision(AABBpair(1:3),AABBpair(4:6),AABBpair(7:9),AABBpair(10:end));

                    if isColliding
                        collisionPairs = [collisionPairs;facesCellIndexObject1,facesCellIndexObject2];
                        AABBpairs = [AABBpairs;AABBpair];
                    end
                    
                end
            end

            for pair = 1:size(collisionPairs,1)
                [JcOut, phiOut, cInfoOut] = findMeshContacts( obj, mesh3D, collisionPairs(pair,1), collisionPairs(pair,2), p,AABBpairs(pair,7:end));
                %meshes(m1), meshes(m2), mdofs{m1}, mdofs{m2} );                    
                if ( isempty(phiOut) )
                    continue;
                end
                Jc = [ Jc; JcOut ];
                phi = [ phi; phiOut ];
                cInfo = [ cInfo, cInfoOut ];

                if ( obj.quality )
                    [JcOut, phiOut, cInfoOut] = findMeshContacts( obj, mesh3D, collisionPairs(pair,2), collisionPairs(pair,1), p, AABBpairs(pair,1:6));
                    %meshes(m1), meshes(m2), mdofs{m1}, mdofs{m2} );                    
                    if ( isempty(phiOut) )
                        continue;
                    end
                    Jc = [ Jc; JcOut ];
                    phi = [ phi; phiOut ];
                    cInfo = [ cInfo, cInfoOut ];
                end
            end
        end
        
        function render(~,~) 
        end
        
        function [ Jc, phi, cInfo ] = findMeshContacts( obj, mesh3D, facesCellIndexObject1, facesCellIndexObject2, p, AABB2_in)
            % Compute contact Jacobians for points on the boundary of mesh1
            % that fall within the boundary of mesh2.
            
            % TODO: looks dangerous to concatenate edges... but we will
            % only do this to get the boundary points!
            % bes1edges = mesh3D.facesSets{facesCellIndexObject1};
            bes1Vert = 1:mesh3D.mergeMeshVerticeCount(facesCellIndexObject1);
            bes2edges = mesh3D.facesSets{facesCellIndexObject2};
            bes2Vert = 1:mesh3D.mergeMeshVerticeCount(facesCellIndexObject2);
            
            adjustedBes1 = bes1Vert+ mesh3D.facesSetsOffset(facesCellIndexObject1);
            adjustedBes2 = bes2Vert+ mesh3D.facesSetsOffset(facesCellIndexObject2);
            
            s1BoundaryPoints     = [ p(adjustedBes1*3-2), p(adjustedBes1*3-1), p(adjustedBes1*3) ];
            s2BoundaryPoints     = [ p(adjustedBes2*3-2), p(adjustedBes2*3-1), p(adjustedBes2*3) ];

            %AABB for vertices too so we can eliminate a bunch of them
            AABB2 = reshape(AABB2_in,3,2);
            s1InAABB2 = ...
                s1BoundaryPoints(:,1) >= AABB2(1,1) & ...
                s1BoundaryPoints(:,1) <= AABB2(1,2) & ...
                s1BoundaryPoints(:,2) >= AABB2(2,1) & ...
                s1BoundaryPoints(:,2) <= AABB2(2,2) & ...
                s1BoundaryPoints(:,3) >= AABB2(3,1) & ...
                s1BoundaryPoints(:,3) <= AABB2(3,2) ;
            in = s1InAABB2;
            if sum(in) <= 0
                Jc = [];
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            adjustedBes1= adjustedBes1(in);
            % Find phi as minimum distance, and choose opposite part of 
            % contact as that closest point on the edge of the other 
            % geometry.  This will violate equal and opposite friction 
            % forces, but is but probably OK if there isn't a lot of
            % interpenetration. 
            m1ContactPoints = s1BoundaryPoints( in, : ); % contact points on bondary of mesh 1

            [S,I,C,N] = signed_distance(m1ContactPoints, s2BoundaryPoints, double(bes2edges));

            % these interpenetraiton depths must be negative for
            % baumgarte to work!
            in = S < 0 & ~isnan(N(:,1)) & ~isnan(N(:,2)) & ~isnan(N(:,3));
            
            m1ContactPoints = m1ContactPoints(in,:);
            
            if sum(in) <= 0
                Jc = [];
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            phi = S(in);
            phi(phi > 0) = 0;
            
            normals = N(in,:);
            tangents = zeros(size(normals));
            tangents2 = zeros(size(normals));
            for i = 1: size(normals,1)
                tangent = null(normals(i,:));
                tangents2(i,:) = tangent(:,2)';
                tangents(i,:) = tangent(:,1)';
            end
            numContacts = sum(in);
            
            % To build Jacobians, note that m1ContactPoints come directly 
            % from DOF indices via their index, while m1ContactPoints come
            % from the nodes of boundary edges.

            % Here, indices are of edge corresponding to mesh1 boundary edges
            % and recall the point was taken as the first index of the edge!
            m1ContactVertIndex = adjustedBes1(in);
            C = C(in,:);
            I = I(in);

            row = 1;
            cInfo = contactInfo3D.empty;
            
            ii = zeros(numContacts*36,1);
            jj =  zeros(numContacts*36,1);
            vals =  zeros(numContacts*36,1);
            for i = 1:numContacts
                verticesOnfaceHit = bes2edges(I(i),:);
                verticesAdjustedIds = adjustedBes2(verticesOnfaceHit);
                points = [p(verticesAdjustedIds*3-2),p(verticesAdjustedIds*3-1),p(verticesAdjustedIds*3)];
                dist = zeros(size(points));
                for j = 1:3
                    dist(:,j) = points(j,:) - C(i,:);
                end
                magnitude = sqrt(sum(dist.^2,2));
                
                alphas = 1 - magnitude/sum(magnitude);
                
                %velocity
                ii(i*36-35:i*36-24) = row;
                jj(i*36-35:i*36-24) = [m1ContactVertIndex(i)*3;m1ContactVertIndex(i)*3-1;m1ContactVertIndex(i)*3-2;
                                      verticesAdjustedIds(1)*3;verticesAdjustedIds(1)*3-1;verticesAdjustedIds(1)*3-2;
                                      verticesAdjustedIds(2)*3;verticesAdjustedIds(2)*3-1;verticesAdjustedIds(2)*3-2;
                                      verticesAdjustedIds(3)*3;verticesAdjustedIds(3)*3-1;verticesAdjustedIds(3)*3-2];
                vals(i*36-35:i*36-24) = [normals(i,3);normals(i,2);normals(i,1);
                                        -normals(i,3)' *alphas(1);-normals(i,2)' *alphas(1);-normals(i,1)' *alphas(1);
                                        -normals(i,3)' *alphas(2);-normals(i,2)' *alphas(2);-normals(i,1)' *alphas(2);
                                        -normals(i,3)' *alphas(3);-normals(i,2)' *alphas(3);-normals(i,1)' *alphas(3);];
                row = row + 1;
                
                %friction
                ii(i*36-23:i*36-12) = row;
                jj(i*36-23:i*36-12) = [m1ContactVertIndex(i)*3;m1ContactVertIndex(i)*3-1;m1ContactVertIndex(i)*3-2;
                                      verticesAdjustedIds(1)*3;verticesAdjustedIds(1)*3-1;verticesAdjustedIds(1)*3-2;
                                      verticesAdjustedIds(2)*3;verticesAdjustedIds(2)*3-1;verticesAdjustedIds(2)*3-2;
                                      verticesAdjustedIds(3)*3;verticesAdjustedIds(3)*3-1;verticesAdjustedIds(3)*3-2];
                vals(i*36-23:i*36-12) = [tangents(i,3);tangents(i,2);tangents(i,1);
                                        -tangents(i,3)' *alphas(1);-tangents(i,2)' *alphas(1);-tangents(i,1)' *alphas(1);
                                        -tangents(i,3)' *alphas(2);-tangents(i,2)' *alphas(2);-tangents(i,1)' *alphas(2);
                                        -tangents(i,3)' *alphas(3);-tangents(i,2)' *alphas(3);-tangents(i,1)' *alphas(3);];
                row = row + 1;
                
                                %friction
                ii(i*36-11:i*36) = row;
                jj(i*36-11:i*36) = [m1ContactVertIndex(i)*3;m1ContactVertIndex(i)*3-1;m1ContactVertIndex(i)*3-2;
                                      verticesAdjustedIds(1)*3;verticesAdjustedIds(1)*3-1;verticesAdjustedIds(1)*3-2;
                                      verticesAdjustedIds(2)*3;verticesAdjustedIds(2)*3-1;verticesAdjustedIds(2)*3-2;
                                      verticesAdjustedIds(3)*3;verticesAdjustedIds(3)*3-1;verticesAdjustedIds(3)*3-2];
                vals(i*36-11:i*36) = [tangents2(i,3);tangents2(i,2);tangents2(i,1);
                                        -tangents2(i,3)' *alphas(1);-tangents2(i,2)' *alphas(1);-tangents2(i,1)' *alphas(1);
                                        -tangents2(i,3)' *alphas(2);-tangents2(i,2)' *alphas(2);-tangents2(i,1)' *alphas(2);
                                        -tangents2(i,3)' *alphas(3);-tangents2(i,2)' *alphas(3);-tangents2(i,1)' *alphas(3);];
                row = row + 1;
                
                % could sort these... but these will likewise always
                % be found in the same order
                indices = [m1ContactVertIndex(i), verticesAdjustedIds];
                cInfo(i) = contactInfo3D([(m1ContactPoints(i,:) + C(i,:))/2], ...
                    normals(i,:), ...
                    tangents(i,:), ...
                    obj.FrictionCoefficient, ...
                    indices, ...
                    obj.ID, ...
                    tangents2(i,:),...
                    phi(i), ...
                    obj.xpbdContactCompliance);
            end   
            Jc = sparse(ii,jj,vals,numContacts*3, 3*mesh3D.N);
        end
    end
end