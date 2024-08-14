classdef NoRigidificator3D < Rigidificator3D
    %EDotMexRigidificator Uses EDot values to tell if a triangle should be
    %made rigid or not. This version uses mex to speedup the process
    properties
        ScaleByMaxEdgeLength
    end
    
    methods
        function obj = NoRigidificator3D()
            %Mexed rigidificator, only uses the elastification call to do
            %both rigidification and elastification in an efficient BFS.
            obj@Rigidificator3D();
        end
        
        function checkForElastification( obj, mesh3D, cache, frame, h, settings )
            
        end
    end
end

