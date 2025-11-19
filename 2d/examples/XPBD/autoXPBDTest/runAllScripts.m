function runAllScripts(path_out, callbackTest, errorStop)
    runOnGPU = false;
    for errorStopTest = errorStop
        callbackTest(path_out, errorStopTest, runOnGPU);
    end
end