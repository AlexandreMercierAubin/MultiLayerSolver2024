classdef Simulation3DSettings < handle
    %Simulation3DSettings stores many simulation parameters
    properties
        %rendering options
        MakeVideo = 0;
        WriteOBJs = 0;
        VideoOut = strcat("out",filesep);
        ExportNormalsOBJs = 0;
        OBJDir = './objs/';
        FramesToRecord = -1;
        SceneName = 'DefaultSceneName';
        PlotSkip = 0;  % plot every frame by default
        MaximizeWindow = 0;
        CamPadding = [0.5 0.5 1 0]; % left, right, down, up
        InitialWindowPosition = [0, 0, 600, 500];
        projection = "perspective";
        campos= [1,0,0];
        camtarget= [0,0,0];
        camfov = 25;
        backgroundColor = 'w';
        renderer = 'zbuffer'; % can also be 'opengl' or 'painters'
        camLightPosition = 'right'; %right, left, or headlight
        camLightStyle = 'infinite'; %infinite or local
        elementLineThickness = 1;
        elementLineColor = 'none';
        
        %simulation options
        RunRightAway = 1; 
        FocusOnMeshNode = [];
        WarmStartEnabled = true;
        PGSiterations = 100;
        
        %visual debugging options
        DrawRigidFrames = 0;
        DrawForces = 0;
        DrawVelocities = 0;
        DrawLagrangeMultipliers = 0;
        DrawEDots = 0;
        DrawTimings = 0;
        DrawLegend = 0;
        DrawEigensOfE = 0;
        DrawEdges = true;
        DrawContact = 0;
        DrawLambdas = 0;
        DrawLambdasScale = 1
        DrawRigidDOFs = 0;
        PlotEDotHist = 0;
        PlotEigs = 0;
        PlotDihedralRateHist = 0;
        PlotEdotVsCurvatureHists = 0;
        PlotPolarDecomposition = 0;
        plotRigidification = 1;
        PlotSpyA = 0
        spyAfig = 0;
        PlotMomentum = 0
        PrintTimings = 0;
        printFrameNumber = false;
        plotFrameTimeHistogram = false;
        residualComparisonLabels = [];
        residualHistBinSize = 0;
        skipPreProjectionResidual = true;
        overwriteComparisonPositionFromX = 0;
        
        %rigidification options
        ElasticContacts = true;
        StrainLimitingEnabled = false;
        ElastificationEnabled = true;
        RigidificationEnabled = true;
        FirstFrameRigidification = false;
        PCGiterations = 1;
        Shading = true;
        quicksolveSimulation = 0;
        runQuickSolve = true;
        
        addShellNormalDeformation = 1;
        addBendingEnergy = 0;
        
        %timing data options
        RecordFramePositionInTD = false;
        
        %cache options
        recomputeCacheAinv = false;
        sameComparisonAinv = false;

        StopAtImprovementPercent = 0;
        
    end
    methods
        function obj = Simulation3DSettings()
            %tries to read a custom simulation settings so users can all
            %have their independant configs.
            if isfile("3d/customSettings3D.m")
                customSettings3D(obj)
            end
        end
    end
end

