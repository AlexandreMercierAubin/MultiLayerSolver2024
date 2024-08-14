classdef Rigid1D
    %RIGID1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        r
        com
        mass
        indices
        vertex2ri
    end
    
    methods
        function obj = Rigid1D(pos,masses,indices)
            obj.mass = sum(masses);
            obj.com = sum(pos.*masses)/obj.mass;
            obj.r = pos - obj.com;
            obj.indices = indices;
            numArray = 1:numel(obj.indices);
            obj.vertex2ri = dictionary(obj.indices,numArray');
        end

        function pos = updateCom(obj, deltaCom)
            obj.com = obj.com+deltaCom;
            pos = obj.com + obj.r;
        end

        function pos = queryVertices(obj, vertexID)
            rID = obj.vertex2ri(vertexID);
            pos = obj.com + obj.r(rID);
        end
    end
end

