
path_out = "./out/auto/";
errorTests = [90];
runAllScripts(path_out,@auto_inflateBox_xpbdLayer,errorTests);
TDhistograms("inflateBox_layer",path_out+"inflateBox_layer/",errorTests,0.0025);