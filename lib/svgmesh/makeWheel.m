
close all
clear;

numTeeth1 = 21;
PA = 5;            % pressure angle (degrees)
Pitch = 3;       % teeth pitch (?)
lRes = 0.7;         % linear resolution for discretization
ACDFac = 1;         % Addendum Circle Diameter factor
DCDFac1 = 1;      % Dedendum Circle Diameter factor (default = 2)

[g1x, g1y, PCD1, ACD1] = MakeGear1( numTeeth1, PA, Pitch, lRes, ACDFac, DCDFac1 );

figure(1);
clf
axis equal;
grid on; hold on;
%plot(g1x, g1y + PCD1/2,'.r');
plot( g1x, g1y ,'-r' );

% Remove the duplicate geometry point before creating geometry
g1x = g1x(1:end-1);
g1y = g1y(1:end-1);

createGearELE( "2d/data/wheel", g1x, g1y, 5, lRes );

return

function createGearELE( filename, g1x, g1y, radius, lRes )
    N = length( g1x );

    N1 = ceil( 2*pi*radius / lRes );
    theta = linspace( 2*pi, 0, N1 );
    theta = theta(1:end-1);
    xdata1 = radius*cos(theta);
    ydata1 = radius*sin(theta);
    N1 = numel(theta);

    % combine the edges of the gear and hole, and make edges
    points = [ g1x, xdata1; g1y, ydata1 ]';
    edges2 = [ 1:N; 2:N,1 ]; 
    edges1 = [ 1:N1; 2:N1,1 ] + [N;N];
    edges = [ edges2, edges1 ]';

    writePOLY_triangle( filename + ".poly", points, edges, []);

    system( pwd + "/lib/svgmesh/triangle.exe  -p -q30 -a0.5 " + filename + ".poly" );
    [TV,I] = readNODE( filename + ".1.node" );
    [TF,A] = readELE( filename + ".1.ele");

    figure
    if numel(A) > 0
        colours = [ 0.5 0.5 0.5; 1 1 0.5; 0.5 1 1; 1 0.5 1; 1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5 ];
        faceColor = colours(mod(A,7)+1,:);
        patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', 'flat', 'FaceVertexCData', faceColor, 'FaceAlpha', .5, 'EdgeAlpha', .9 );
    else
        patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', 1 );

    end
    axis equal
    axis off
    
end

