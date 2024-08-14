
path_out = "./out/auto/";
% errorTests = [30,60,90];
% errorTests = [70,80,90];
errorTests = [90];
% 
% runAllScripts(path_out,@auto_cantilever_xpbdLayer, errorTests);
% TDhistograms("cantilever_xpbdLayer",path_out+"cantilever_xpbdLayer/",errorTests);
% runAllScripts(path_out,@auto_spider_xpbdLayer,errorTests);
% TDhistograms("spider_xpbdLayer",path_out+"spider_xpbdLayer/",errorTests);
% runAllScripts(path_out,@auto_IShape_xpbdLayer,errorTests);
% TDhistograms("IShapeGravity_xpbdLayer",path_out+"IShapeGravity_xpbdLayer/",errorTests);
% runAllScripts(path_out,@auto_crossPinned,errorTests);
% TDhistograms("crossPinned",path_out+"crossPinned/",errorTests);

% runAllScripts(path_out,@auto_inflateBox_xpbdLayer,errorTests);
% TDhistograms("inflateBox_layer",path_out+"inflateBox_layer/",errorTests);

runAllScripts(path_out,@auto_ground_wheel,errorTests);
TDhistograms("groundWheel",path_out+"groundWheel/",errorTests,0.005);

% runAllScripts(path_out,@auto_inflateIShapeCenter,errorTests);
% TDhistograms("inflateIShape_layer",path_out+"inflateIShape_layer/",errorTests,0.0005);
% 
% runAllScripts(path_out,@auto_spinBox_xpbdLayer,errorTests);
% TDhistograms("spinBox_layer",path_out+"spinBox_layer/",errorTests,0.02);

% runAllScripts(path_out,@auto_octopus_XPBD, errorTests);
% TDhistograms("auto_octopus_xpbd","./out/auto/auto_octopus_xpbd/",errorTests,0.01,true);
% 
% run("cantilever_xpbdLayerOrder.m");
% run("cantilever_xpbdLayerCycles.m");
% run("cantilever_xpbdLayerNumberOfLayers.m");
% run("cantilever_xpbdLayerNumberIterations.m");
% run("cantilever_xpbdLayerSolversComp.m");
% run("primitivesPachinkoFull.m");
% TDhistograms("PachinkoPill","./",[""],0.06);
% run("octopus_XPBD.m");
% TDhistograms("octopusPlatform","./",[""],0.01);


function runAllScripts(path_out, callbackTest, errorStop)
    runOnGPU = false;
    for errorStopTest = errorStop
        callbackTest(path_out, errorStopTest, runOnGPU);
    end
end