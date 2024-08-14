classdef PlatformContactDetector < ContactFinder
    %PlatformContactDetector collision detector that allows the collision
    %detection between meshes and a rectangular prism platform
    
    properties
        HalfWidth
        HalfHeight
        HalfDepth
        DirectionX
        DirectionY
        DirectionZ
        Direction
        Center
        V %vertices
        F %faces
        R
        TangentMapping = [2,3;
                          1,3;
                          1,2];
        Color = 'red';
    end
    
    methods
        function obj = PlatformContactDetector( degrees, scale, position, frictionCoefficient )
            obj.HalfWidth = scale(1)/2;
            obj.HalfDepth = scale(2)/2;
            obj.HalfHeight = scale(3)/2;

            obj.Center = position;
            
            obj.R = RotationMatrixDegrees(degrees);
            
            obj.DirectionX = obj.R(1,:)./norm(obj.R(1,:),2);
            obj.DirectionY = obj.R(2,:)./norm(obj.R(2,:),2);
            obj.DirectionZ = obj.R(3,:)./norm(obj.R(3,:),2);
            obj.Direction = [obj.DirectionX;obj.DirectionY;obj.DirectionZ];
            
            if nargin >= 8
                obj.FrictionCoefficient = frictionCoefficient;
            end
            
            obj.V(1,1:3) = obj.Center - obj.DirectionX * obj.HalfWidth - obj.DirectionY * obj.HalfDepth - obj.DirectionZ * obj.HalfHeight;
            obj.V(2,1:3) = obj.Center - obj.DirectionX * obj.HalfWidth - obj.DirectionY * obj.HalfDepth + obj.DirectionZ * obj.HalfHeight;
            obj.V(3,1:3) = obj.Center + obj.DirectionX * obj.HalfWidth - obj.DirectionY * obj.HalfDepth - obj.DirectionZ * obj.HalfHeight;
            obj.V(4,1:3) = obj.Center + obj.DirectionX * obj.HalfWidth - obj.DirectionY * obj.HalfDepth + obj.DirectionZ * obj.HalfHeight;
            obj.V(5,1:3) = obj.Center - obj.DirectionX * obj.HalfWidth + obj.DirectionY * obj.HalfDepth - obj.DirectionZ * obj.HalfHeight;
            obj.V(6,1:3) = obj.Center - obj.DirectionX * obj.HalfWidth + obj.DirectionY * obj.HalfDepth + obj.DirectionZ * obj.HalfHeight;
            obj.V(7,1:3) = obj.Center + obj.DirectionX * obj.HalfWidth + obj.DirectionY * obj.HalfDepth - obj.DirectionZ * obj.HalfHeight;
            obj.V(8,1:3) = obj.Center + obj.DirectionX * obj.HalfWidth + obj.DirectionY * obj.HalfDepth + obj.DirectionZ * obj.HalfHeight;
            
            obj.F = [1,6,5;
                     1,2,6;
                     2,8,6;
                     2,4,8;
                     4,7,8;
                     4,3,7;
                     3,5,7;
                     3,1,5;
                     6,8,7;
                     6,7,5;
                     2,3,4;
                     2,3,1]; 
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time, p)
            ps = vertcat(p);
            vertices = formatPositions3D(ps);
            directions = vertices-obj.Center;
            
            % indices of the points within the box
            dotX = dot(directions, repmat(obj.DirectionX,size(directions,1),1),2);
            absDotX = abs(dotX);
            dotY = dot(directions, repmat(obj.DirectionY,size(directions,1),1),2);
            absDotY = abs(dotY);
            dotZ = dot(directions, repmat(obj.DirectionZ,size(directions,1),1),2);
            absDotZ = abs(dotZ);
            idx = absDotX <= obj.HalfWidth & ...
                    absDotY <= obj.HalfDepth & ...
                    absDotZ <= obj.HalfHeight;
            
            if ( all(idx == false) )
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            axesDistances = [obj.HalfWidth-absDotX(idx), obj.HalfDepth - absDotY(idx), obj.HalfHeight - absDotZ(idx)];
            [phiDist,directionCase] = min(axesDistances,[],2); %rowwise min
            phi = -phiDist;
            %directionCase is the orientation axis which is the closest
            %X,Y,Z
            contactDots = [dotX(idx),dotY(idx),dotZ(idx)];
            contactSignsVector=sign(contactDots(sub2ind(size(contactDots),1:size(contactDots,1),directionCase')));
            normal = contactSignsVector'.*obj.Direction(directionCase,:);

            %Tangents can just be the box orientations that are farther regardless of signs
            tangent = obj.Direction(obj.TangentMapping(directionCase,1),:);
            tangent2 = obj.Direction(obj.TangentMapping(directionCase,2),:);
            
            indices = find(idx);

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
                cInfo(i) = contactInfo3D( vertices(indices(i),:), normal(i,:), tangent(i,:), obj.FrictionCoefficient, indices(i), obj.ID, tangent2(i,:),phi(i), obj.xpbdContactCompliance);
            end             
        end
        
        function render( obj, ~ )
            if ( obj.plotHandle ~= 0 )
                return
            end
            
            hold on;
            obj.plotHandle = patch('Faces',obj.F,'Vertices',obj.V,'FaceColor',obj.Color);
        end
        
        function [V,F] = getObjPositionFaces(obj, time)
            F = obj.F;
            V = obj.V;
        end
    end
end

