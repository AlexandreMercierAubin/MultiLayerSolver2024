function makeZBar()
    lres = 0.2;
    maxArea = 0.03;    
    filename = "2d/data/ZBar";
    
    points = [-1.5,-1.5;
         -1,  -1.5;
         -1,  -1;
          2,  -1;
          2,   -0.5;
          1.5, -0.5;
          1.5, 0;
          1,   0;
          1,   -0.5;
          -2,  -0.5;
          -2, -1;
          -1.5, -1]; 
    ms = 15; 
    figure(1);
    clf
    
    plot( points(:,1),points(:,2), 'r.', 'Markersize', ms );
    title( "lres = " + lres );

    % Save this to a poly file, and triangulate!
    N = size(points,1);
    edges = [ 1:N; 2:N,1 ]'; 

    writePOLY_triangle( filename + ".poly", points, edges, []);    

    system( pwd + "/lib/svgmesh/triangle.exe  -p -q30 -a"+ maxArea + " " + filename + ".poly" );
    [TV,I] = readNODE( filename + ".1.node" );
    [TF,A] = readELE( filename + ".1.ele");

    figure
    patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .9 );
    axis equal
    axis off
end   