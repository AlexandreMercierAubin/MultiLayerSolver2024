classdef AnimatedPart < handle
    %ANIMATEDPART Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        VerticesRest
        Faces
        Positions    %n by 3 position of the reference frame
        Orientations %n by 4 quaternion orientation of the part
    end
    
    methods
        function obj = AnimatedPart(folder, entry, positions, orientations, scale)
            filename = [folder,entry{1},num2str(entry{2},'%04.f'),'.obj'];
            [obj.VerticesRest,obj.Faces,UV,TF,N,NF] = readOBJ(filename);
            obj.VerticesRest = (diag(scale)*obj.VerticesRest')';
            obj.Positions = positions;
            obj.Orientations = orientations;
        end

        function vertices = getVerticesAtFrame(obj, frame, scaleRatio)
            refFrame = frame;
            if frame > size(obj.Orientations,1)
                refFrame = size(obj.Orientations,1);
            end
            
            R = quat2mat(obj.Orientations(refFrame,:));

            rotatedVerts = (R*scaleRatio*obj.VerticesRest')';
            vertices = rotatedVerts + obj.Positions(refFrame,:);

        end
        
    end
end

