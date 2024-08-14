classdef SimulationSettings < handle
    % Settings for a simulation, including viewing, drawing options, and 
    % data collection options, and others.  See 'doc SimulationSettings' 
    % for more information on individual fields.
    properties
        %rendering options
        MakeVideo = 0;
        VideoOut = strcat("out",filesep);
        ExportNormalsOBJs = 0;
        FramesToRecord = -1;
        SceneName = 'DefaultSceneName'; % scene name is important as it will be used for both video output filename and caching.
        PlotSkip = 0;  % plot every frame by default
        MaximizeWindow = 0;
        CamPadding = [0.5 0.5 1 0]; % left, right, down, up
        InitialWindowPosition = [500, 200, 600, 500]; %initial position of the window x,y on the screen, width, height
        projection = "perspective";
        Shading = true;
        campos= [0,0,0];                             %position of the camera in the scene, parameter x,y are the important ones in 2D
        camtarget= [0,0,0];
        camfov = 25;
        MakeWindowSquareBounds = true;
        elementLineThickness = 1;
        elementLineColor = [0,0,0];
        boundaryColor = 'black';
        boundaryThickness = 1.5;
        elementLineAlpha = 0.05;
        elementFaceAlpha = 0.5; %TODO: maybe move this to the material properties with the rest of the color so we have per element alphas
        backgroundColor = 'w';

        %simulation options
        RunRightAway = 1;                           %start the scene immediately or press space to start
        FocusOnMeshNode = [];                       %allows the camera to track a vertex of a mesh
        WarmStartEnabled = true;                    %warm starts the contact solve with previous lambdas
        PGSiterations = 50; % Projected Gauss-Seidel iterations for contact solve.
        StrainLimitingEnabled = false;
        overwriteComparisonPositionFrom1 = false;
        StopAtImprovementPercent = 0;
        
        %visual debugging options
        DrawGraphColors = 0;
        DrawContactFrames = 0;
        DrawRigidFrames = 0;
        DrawForces = 0;
        DrawVelocities = 0;
        DrawLagrangeMultipliers = 0;
        DrawEDots = 0;
        DrawApproxEDots = 0;
        DrawTimings = 0;
        DrawLegend = 0;
        DrawEigensOfE = 0;
        DrawEdges = true;
        DrawContact = 0;
        DrawLambdas = 0;
        DrawLambdasScale = 50;
        DrawRigidDOFs = 0;
        DrawApproxDv = false;
        DrawDv = false;
        DrawElasticForces = false;
        PlotEDotHist = 0;
        PlotRateWork = 0;
        PlotApproxRateWork = 0;
        PlotPhiHist = 0;
        plotList = [];
        PlotPolarDecomposition = 0;
        plotTriImplicitRigidificationElastification = 0;
        plotRigidification = 1;
        PlotContactForceProfile = 0;
        PlotSpyA = 0;
        PopUpPlot = false;
        spyAfig
        PlotMomentum = 0;
        PrintTimings = 0;
        plotLayers = 0;
        plotOnlyComparison = -1;
        plotLayerSwitch = false;
        plotResidual = false;
        plotFrameTimeHistogram = false;
        residualComparisonLabels = ["Elastic","Multi-Layers"];
        plotResidualRelativeComparison = false;
        plotResidualFrames = inf;%min between the number of recorded frames and this for the plots.
        residualRange = [-inf,inf];
        residualHistBinSize = -1; %uses the first bin width unless specified
        skipPreProjectionResidual = true;
        specialResidualOption = false;
        StopCompareToID = [1];
        
        %rigidification options
        ElasticContacts = true;                         %quick parameter to disable meshmesh contacts
        ElastificationEnabled = true;                   
        RigidificationEnabled = true;
        FirstFrameRigidification = true;
        PCGiterations = 1; % Preconditioned conjugate gradient iterations for quick solve. 
        quicksolveSimulation = 0; %use solely to visualize what the quicksolve outputs
        warmStartFromQuickSolve = false; %use the lambdas from the quicksolve contact handling to warmstart the pgs
        plotElasticIntegrationTime = false;
        
        %timing data options
        RecordFramePositionInTD = false;            %keep an array of point position in the timing data
        
        %cache options
        recomputeCacheAinv = false; % forces recompute of the Ainv blocks and the preconditioner
        sameComparisonAinv = false;
    end
    methods
        function obj = SimulationSettings(customSettingsFunc)
            %tries to read a custom simulation settings so users can all
            %have their independant configs.
            if nargin > 0
                obj = customSettingsFunc(obj);
            end
        end
    end
end

