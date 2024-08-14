classdef AnimatedObjectContactDetector < ContactFinder
    %ObjectContactDetector collision detector that allows the collision
    %detection between meshes and an imported obj mesh
    
    properties
        parts
        skipdetection = false;
        interpolateOnFrames = 1;
        plotScale = 1;
        edgeColor = 'none';
    end
    
    methods
        function obj = AnimatedObjectContactDetector(filename, partsFolder, degrees, scale, position, frictionCoefficient)
            [positions,orientations,entries] = readMOT(filename);
            obj.parts = [];

            angles = deg2rad(degrees);
            quatAngles = angle2quat(angles(1),angles(2),angles(3));
            R = quat2mat(quatAngles);

            for i = 1:numel(entries)
                scaledPositions = scale.*positions(:,:,i);
                transformedPositions = (R*scaledPositions')' + position;
                rotatedQuats = quatmultiply(quatAngles,orientations(:,:,i));
                obj.parts = [obj.parts, AnimatedPart(partsFolder,entries{i},transformedPositions , rotatedQuats, scale)];
            end
            
            obj.FrictionCoefficient = frictionCoefficient;
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time, p)
            ps = vertcat(p);
            xs = ps(1:3:end);
            ys = ps(2:3:end);
            zs = ps(3:3:end);
            
            [V,F] = obj.getObjPositionFaces(time,1);
            
            % AABB
            [minV,maxV] = bounds(V,1);
            AABB = [minV',maxV'];
            s1InAABB2 = ...
                xs >= AABB(1,1) & ...
                xs <= AABB(1,2) & ...
                ys >= AABB(2,1) & ...
                ys <= AABB(2,2) & ...
                zs >= AABB(3,1) & ...
                zs <= AABB(3,2) ;
            idx = s1InAABB2;

            if ( sum(idx) == 0  || obj.skipdetection)
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            points = [xs(idx),ys(idx),zs(idx)];
            
           [S,I,C,N] = signed_distance(points, V, F);

            % these interpenetraiton depths must be negative for
            % baumgarte to work!
            in = S < 0 & ~isnan(N(:,1)) & ~isnan(N(:,2)) & ~isnan(N(:,3));
            
            points = points(in,:);
            
            if sum(in) <= 0
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            phi = S(in);
            phi(phi > 0) = 0;
            
            normal = N(in,:);
            tangent = zeros(size(normal));
            tangent2 = zeros(size(normal));
            for i = 1: size(normal,1)
                tangents = null(normal(i,:));
                tangent(i,:) = tangents(:,1)';
                tangent2(i,:) = tangents(:,2)';
            end
            
            indices = find(idx);
            indices = indices(in);

            %add friction
            rown  = (1:3:3*numel(phi))';
            rowt  = (2:3:3*numel(phi))';
            rowt2  = (3:3:3*numel(phi))';
            colx = (indices*3-2);
            coly = (indices*3-1);
            colz = (indices*3);
            J = sparse( ...
                [ rown;         rown;      rown;       rowt;         rowt;     rowt;   rowt2;         rowt2;     rowt2], ...
                [ colx;         coly;      colz;       colx;         coly;     colz;   colx;         coly;     colz;], ...
                [ normal(:,1);  normal(:,2); normal(:,3); tangent(:,1); tangent(:,2) ; tangent(:,3) ; tangent2(:,1); tangent2(:,2) ; tangent2(:,3)], 3*numel(phi), numel(ps) );

            cInfo = contactInfo3D.empty; %cell( numel(phi), 1 );
            for i = 1:numel(phi) 
                cInfo(i) = contactInfo3D( [xs(indices(i)), ys(indices(i)), zs(indices(i))], normal(i,:), tangent(i,:), obj.FrictionCoefficient, indices(i), obj.ID, tangent2(i,:),phi(i));
            end             
        end
        
        function render( obj, frame )
            [V, F] = obj.getObjPositionFaces(frame, obj.plotScale);
            
            if ( obj.plotHandle ~= 0 )
                obj.plotHandle.Vertices = V;
                return
            end

            hold on;
            obj.plotHandle = patch('Faces',F,'Vertices',V,'FaceColor','red', 'EdgeColor', obj.edgeColor);
        end
        
        function [V,F] = getObjPositionFaces(obj, frame, scaleRatio)
            V = [];
            time = floor(frame/obj.interpolateOnFrames)+1;
            modValue = mod(frame,obj.interpolateOnFrames);
            alpha = modValue / obj.interpolateOnFrames;
            for i = 1:numel(obj.parts)
                if frame == 1
                    V = [V;obj.parts(i).getVerticesAtFrame(time,scaleRatio)];
                else
                    V0 = obj.parts(i).getVerticesAtFrame(time,scaleRatio);
                    V1 = obj.parts(i).getVerticesAtFrame(time+1,scaleRatio);
                    
                    Vobject = alpha*V1 + (1-alpha)*V0;
                    V = [V;Vobject];
                end
            end

            Fcells = {obj.parts.Faces};
            F = Fcells{1};
            for i = 2:numel(obj.parts)
                maxF = max(F(:));
                F = [F;Fcells{i}+maxF];
            end
            
        end
    end
end

