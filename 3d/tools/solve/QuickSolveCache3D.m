classdef QuickSolveCache3D < handle
    % QUICKSOLVECACHE Cache object that holds values that are computed in
    % the quickSolve step and then recycled by the integrator
    
    properties
        ApproximatedDeltaV
        recomputeGravityDv = 1;
        gravityDv
        cInfo = {}
        contactIDs = []   % maintained here and in the cInfo :(  Here is more useful
        
        prevCInfo = {}
        prevLambdas = []
        prevContactIDs = []  % maintained here and in the cInfo :(  Here is more useful
        
        AInvBlocks = []   % computed on the first pass through... could be saved/loaded
        
        F                 % deformation gradients
        dpsidF            % negative area times energy gradient wrt F for each element
        elasticForces     % elastic forces
        C                 % negative area times energy hessian wrt F
        shellBendingH     % bending energy Hessian of shell
        bendingForces     % bending energy forces of shell
        NormalShellComponent %change in normals for shells
        dihedralApproxRate %for TriMesh3D, bending rate from quicksolve
        dihedralRate       %for TriMesh3D, bending rate
        
        % preconditioner stuff... if we need a closer look later, could put
        % the matrices in here too
        preconditioner    % preconditioner function
        
        warmStartLambdasComputed = false
        wsl
        ncf   
        
        oldp
        oldOldp
        oldv
        oldDv
        oldF
        prevPrevDihedralAngle
        edotnorms
        edotApproxNorms
        Perm
        Lpre
        App
        Apre
        ActiveB % this is caching the B matrix containing active rows as indexing is slow.
        dndp %this is the derivative of the normals with respect to the positions
        layers
        initialResidual
        residualArray
    end
    
    methods
        function obj = QuickSolveCache3D()
        end
        
        function clearWarmStartInfo( obj )
            % clearWarmStartInfo resets the cache as a step has been taken
            % and matching contacts will need to be found again after
            % running collision detection again.  Note that this call coule
            % easily be an update cinfo call as the prevCinfo gets set at
            % the same time as this is called
            obj.warmStartLambdasComputed = false; 
        end
        
        function clear( obj )
            obj.ApproximatedDeltaV = [];
            obj.cInfo = contactInfo.empty;
            obj.prevCInfo = contactInfo.empty;
            obj.prevLambdas = [];
            obj.warmStartLambdasComputed = false;
        end
        
        function [warmStartLambdas, newContactFlag] = findWarmStartLambdaArray( obj )
            if obj.warmStartLambdasComputed
                warmStartLambdas = obj.wsl;
                newContactFlag = obj.ncf;
                return;
            end
            obj.warmStartLambdasComputed = true;
            if ( numel(obj.prevCInfo) == 0 ) 
                warmStartLambdas = zeros([3*numel(obj.cInfo),1]);
                newContactFlag = true([numel(obj.cInfo),1]); 
                obj.wsl = warmStartLambdas;
                obj.ncf = newContactFlag;
                return;
            end

            [ warmStartLambdas, newContactFlag ] = mexFindMatches3D( obj.prevContactIDs, obj.prevLambdas, obj.contactIDs );

            obj.wsl = warmStartLambdas;
            obj.ncf = newContactFlag;
        end
    end
end

