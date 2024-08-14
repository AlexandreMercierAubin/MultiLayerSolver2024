classdef NoRigidificator < Rigidificator
    %EDotMexRigidificator Uses EDot values to tell if a triangle should be
    %made rigid or not. This version uses mex to speedup the process
    properties
        ScaleByMaxEdgeLength
        ForceRigidTriList = [];
        forceRebuild = false;
    end
    
    methods
        function obj = NoRigidificator()
            %EDOTVECTORRIGIDIFICATOR Same as EDotRigidficator but
            %vectorized and thus much faster
            obj@Rigidificator();
        end
        
        function checkForElastification( obj, mesh2D, cache, frame, h, settings )
        end
    end
end

