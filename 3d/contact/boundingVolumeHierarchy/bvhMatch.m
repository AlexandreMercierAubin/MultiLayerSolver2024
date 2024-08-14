classdef bvhMatch
    %BVHMATCH Contains information about two faces that are
            %potentially in contact. This is meant to be created on
            %leaf-leaf collision
    
    properties
        centerDistance,
        faceID1,
        faceID2
    end

    methods
        function obj = bvhMatch(centerDistance, faceID1,faceID2)
            %CHILDBVH Contains information about two faces that are
            %potentially in contact. This is meant to be created on
            %leaf-leaf collision
            %centerDistance: distance between the centers of the bounding
            %volumes
            %faceID1: id of the first face in the collision
            %faceID2: id of the other colliding face
            obj.centerDistance = centerDistance;
            obj.faceID1 = faceID1;
            obj.faceID2 = faceID2;
        end
    end
end

