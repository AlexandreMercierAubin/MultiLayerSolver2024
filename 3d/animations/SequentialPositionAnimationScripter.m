classdef SequentialPositionAnimationScripter < AnimationScripter
    %SequentialPositionAnimationScripter animation scripter that forces
    %nodes to be a specific location at a given frame. The input must be sorted.
    %NOTE that this only works for the last dof of a rigid body

    
    properties
        frameNumbers = [];
        %assumes sorted DOFs
        dofs = {};% n frame to constrain by constrained node at frame
        positions = {}; %position of node index from nodes at frame
        counter = 1;
        dim = 2;
        
    end
    
    methods
        function obj = SequentialPositionAnimationScripter()
            %! NOTE for the mesh to remain stable with no explosion it is
            %recommended to pin all the points that are being moved by the
            %sequential position animation scripter
        end
        
        function obj = scriptMesh(obj, mesh3D, integrator, frame, h)
            return;
        end
        
        function newDV = scriptVelocity(obj, prevP, prevV, deltaV, activeDOFs, activeDOFsIDs, frame, rigidIDbyVert, h)
            newDV = deltaV;
        end
        
        function newDV = scriptVelocityQuicksolve(obj, prevP, prevV, deltaV, frame, h)
            newDV = deltaV;
        end
        
        function [newP,newV] = scriptPositions(obj, prevP, prevV, frame, h)
            % [newP,newV] = scriptPositions(obj, prevP, prevV, frame)
            % prevP: position prior to the script
            % newP: new position of the points of the mesh
            % prevP: previous points position of the mesh
            % frame: frame of the simulation
        
            newP = prevP;
            newV = prevV;
           
            isComparison = (obj.counter>=2 && obj.frameNumbers(obj.counter-1) == frame && obj.counter-1 <= numel(obj.frameNumbers));
            if (obj.counter <= numel(obj.frameNumbers) &&(obj.frameNumbers(obj.counter) == frame) || isComparison)
                frameDOFs = obj.dofs{obj.counter - isComparison};
                framePositions = obj.positions{obj.counter - isComparison};
                newP(frameDOFs) = framePositions;
                newV(frameDOFs) = (newP(frameDOFs) - prevP(frameDOFs))./h;
                if ~isComparison
                    obj.counter = obj.counter + 1;
                end
            end
        end

        function newForce = scriptForceAnimation(obj, force, prevP, prevV,mass, frame, h)
            % not used       
            newForce = force;

        end
        
        function reset(obj)
            obj.counter = 1;
        end
    end
end

