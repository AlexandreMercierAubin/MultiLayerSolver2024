classdef nodeBVH
    %CHILDBVH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        childs
        faces
        faceVerts
        useCubeCollision = true;
    end
    
    methods
        function obj = nodeBVH(faces, faceVerts, childs)
            %CHILDBVH Constructs a BVH node. The top level node (root) should
            %contain different objects to check collision against. Each individual object should have
            %its own tree and be a child of the root. Each level contains a
            %smaller subset of elements.
            obj.childs = childs;
            obj.faces = faces;
            obj.faceVerts = faceVerts;
        end

        function val = isLeaf(obj)
            %validates if this node is at the bottom of the tree (a leaf)
            val = isempty(obj.childs);
        end
        
        function [collides,dist, center, radius] = collisionTest(obj,otherChild, V1, V2)
            %Tests if the faces within two childs collide
            %otherChild: child to check collision against
            %V1: n by 3 current vertex positions of the first child
            %V2: n by 3 current vertex positions of the second child
            fStack = obj.faceVerts(:);
            positions = V1(fStack,:);
            [minPos1,maxPos1] = bounds(positions);

            fStack2 = otherChild.faceVerts;
            positions2 = V2(fStack2,:);
            [minPos2,maxPos2] = bounds(positions2);

            %This uses cubes instead of actual boxes in case a face is axis
            %aligned. This gives more opportunity to detect collisions.
            if obj.useCubeCollision
                [collides,dist, center, radius] = AABBcubeCollision(minPos1, maxPos1, minPos2, maxPos2);
            else
                [collides,dist, center, radius] = AABBboxesCollision(minPos1, maxPos1, minPos2, maxPos2);
            end
        end

        function [collides,dist, center, radius] = collisionTestCCD(obj,otherChild, V1, V2, prevV1, prevV2)
            %Tests if the faces within two childs collide
            %otherChild: child to check collision against
            %V1: n by 3 current vertex positions of the first child
            %V2: n by 3 current vertex positions of the second child
            %prevV1: n by 3 positions of the vertices of the first child at the previous time
            %step
            %prevV2: n by 3 positions of the vertices of the second child at the previous time
            %step

            %TODO: build two prisms from prevV1 to V1 & prevV2 to V2
            %TODO: check prism collision
            %TODO: fill the return values
            collides = false;
            dist = 0;
            center = [0,0,0];
            radius = 0;
        end

        function [collides, matches, toPlot] = intersect(obj,otherChild, V1, V2, isDrawing)
            %This is a recursive function to handle the exploration of the
            %BVH
            %otherChild: child to check collision against
            %V1: n by 3 current vertex positions of the first child
            %V2: n by 3 current vertex positions of the second child
            %isDrawing: plot the bounding volumes?

            %TODO:have an input that decides if we check collision with CCD
            %or with the standard check
            [collides,dist, center, radius] = obj.collisionTest(otherChild, V1, V2);
            matches = [];
            toPlot = [];
            noChildCollision = true; 
            if collides
                if obj.isLeaf() && otherChild.isLeaf() %Two leaves, then stop, we are at the bottom of the tree
                    matches = [bvhMatch(dist,obj.faces, otherChild.faces)];
                    noChildCollision = false;
                    if isDrawing
                        toPlot = [toPlot;center,radius];
                    end
                elseif obj.isLeaf() %one leaf then explore the remaining nodes of the other child
                    for j = 1:numel(otherChild.childs)
                        [childCollision, newMatches, childToPlot] = obj.unevenExplore(otherChild.childs(j), V1, V2, isDrawing);
                        matches = [matches;newMatches];
                        noChildCollision = noChildCollision && ~childCollision;
                        if isDrawing
                            toPlot = [toPlot;childToPlot];
                        end
                    end
                elseif otherChild.isLeaf() %one leaf then explore the remaining nodes of the first child
                    for j = 1:numel(obj.childs)
                        [childCollision, flippedMatches, childToPlot] = otherChild.unevenExplore(obj.childs(j), V2, V1, isDrawing);
                        newMatches = fliplr(flippedMatches);
                        matches = [matches; newMatches];
                        noChildCollision = noChildCollision && ~childCollision;
                        if isDrawing
                            toPlot = [toPlot;childToPlot];
                        end
                    end
                else
                    for i = 1:numel(obj.childs) %no leaves, then explore everything
                        for j = 1:numel(otherChild.childs)
                            [childCollision, newMatches, childToPlot] = obj.childs(i).intersect(otherChild.childs(j), V1, V2, isDrawing);
                            matches = [matches;newMatches];
                            noChildCollision = noChildCollision && ~childCollision;
                            if isDrawing
                                toPlot = [toPlot;childToPlot];
                            end
                        end
                    end
                end
            end
            if noChildCollision && isDrawing
                toPlot = [center,radius];
            end
        end

        function [collides, matches, toPlot] = unevenExplore(obj, otherChild, V1, V2, isDrawing)
            %This is a recursive function to handle the exploration of the
            %BVH when one of the childs is a leaf
            %otherChild: child to check collision against
            %V1: n by 3 current vertex positions of the first child
            %V2: n by 3 current vertex positions of the second child
            %isDrawing: plot the bounding volumes?

            %TODO:have an input that decides if we check collision with CCD
            %or with the standard check
            [collides,dist, center, radius] = obj.collisionTest(otherChild, V1, V2);
            matches = [];
            toPlot = [];
            noChildCollision = true;
            if collides
                if otherChild.isLeaf()% both childs are leaves, we're done exploring
                    matches = [bvhMatch(dist,obj.faces, otherChild.faces)];
                    if isDrawing
                        toPlot = [toPlot;center,radius];
                    end
                    noChildCollision = false;
                else 
                    for j = 1:numel(otherChild.childs) % otherchild is not yet a leaf, then explore its childs
                        [childCollides,childMatches, childToPlot] = obj.unevenExplore(otherChild.childs(j), V1, V2, isDrawing);
                        matches = [matches; childMatches];
                        noChildCollision = noChildCollision && ~childCollides;
                        if isDrawing
                            toPlot = [toPlot;childToPlot];
                        end
                    end
                end
            end
            if noChildCollision && isDrawing
                toPlot = [center,radius];
            end
        end
    end
    methods(Static)
        function child = buildBottomUpNearestHierarchy(V,F, isDrawing)
            %bottomup function, starts with the childs where each box is a face
            %for each level, group boxes to the nearest box in the level
            %keep doing this until you have a single box for an object
            %each object should become a child of the root.
            %make a function to intersect with outside objects like a collision
            %detector
            %this should work as is with both tets and triangles
            
            dim = size(F,2);

            %make each element a node
            [B1,B2] = box_each_element(V,F);
            child = [];
            for i = 1 : size(F,1)
                child = [child,nodeBVH([i],[F(i,:)],[])];
            end

            p = reshape(V', [], 1);
            while numel(child)>1
                childCenter = zeros(numel(child),3);
                for i = 1:numel(child) % computes the distance between bounding volumes
                    Flocal = F(child(i).faces,:);
                    fStack = Flocal(:);
                    positions = V(fStack,:);
%                     [childMin,childMax] = bounds(positions);
%                     childCenter(i,:) = (childMax+childMin)./2;
                    childCenter(i,:) = sum(positions,1)./size(positions,1);
                end

                parents = [];
                while numel(child)>1 %make pairs of closest bounding volumes and group them into a parent
                    lastChildID = numel(child);
                    lastChild = child(lastChildID);
                    dist=sqrt(sum(bsxfun(@minus, childCenter(1:lastChildID-1,:), childCenter(lastChildID,:)).^2,2));
                    nearestChildID = find(dist==min(dist));
                    nearestChildID = nearestChildID(1);
                    nearestChild = child(nearestChildID);
                    
                    faceIDs = [[lastChild.faces];[nearestChild.faces]];
                    faceGroup = [[lastChild.faceVerts];[nearestChild.faceVerts]];
                    childs = [nearestChild,lastChild];
                    newNode = nodeBVH(faceIDs,faceGroup,childs);
                    parents = [parents, newNode];

                    childCenter(lastChildID,:) = [];
                    child(lastChildID) = [];
                    childCenter(nearestChildID,:) = [];
                    child(nearestChildID) = [];
                end
                
                if numel(child) == 1
                    parents = [parents, child];
                end

                child = parents;
            end
        end
    end
end

