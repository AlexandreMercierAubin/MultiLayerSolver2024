
width = 10;
length = 150;
xdata2 = [0,0,length,length];
ydata2 = [0,width,width,0];

shape = polyshape( xdata2, ydata2 );

figure(1);
clf;
subplot(2,1,1);
plot(shape);
hold on;
plot(xdata2,ydata2, '.r' )

axis equal;


filename = "2d/data/long_long_bar";
N = 4;
points = [ xdata2; ydata2 ]';
edges2 = [ 1:N; 2:N,1 ]; 
edges = [ edges2 ]';
holes = [ ];
% points #V by 2
% edges  #E by 2
% holes  #H by 2
writePOLY_triangle( filename + ".poly", points, edges, holes);

system( pwd + "/lib/svgmesh/triangle.exe  -p -q30 -a3 " + filename + ".poly" );
[TV,I] = readNODE( filename + ".1.node" );
[TF,A] = readELE( filename + ".1.ele");

subplot(2,1,2)
if numel(A) > 0
    colours = [ 0.5 0.5 0.5; 1 1 0.5; 0.5 1 1; 1 0.5 1; 1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5 ];
    faceColor = colours(mod(A,7)+1,:);
    patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', 'flat', 'FaceVertexCData', faceColor, 'FaceAlpha', .5, 'EdgeAlpha', .9 );
else
    patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .9 );

end
axis equal
axis off



