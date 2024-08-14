function debugPlotTriangle(p,T)
    V = reshape(p, 2, [])';
    faceColor = repmat([1,0,0],size(T,1),1);
    patch('vertices', V, 'faces', T,'edgecol','black', 'facecol', 'flat','FaceVertexCData', faceColor);
end