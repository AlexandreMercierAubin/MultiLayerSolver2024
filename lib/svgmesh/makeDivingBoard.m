% Make a tube with circular cross section

height = 0.05;
XLen = 1;
YLen = 3;
x=0;
y=0;
z=0;
V = [x, y+YLen/2, z];
V = [V; x, y-YLen/2, z];
V = [V; x+XLen, y+YLen/2, z];
V = [V; x+XLen, y-YLen/2, z];

V = [V; x, y+YLen/2, z+height];
V = [V; x, y-YLen/2, z+height];
V = [V; x+XLen, y+YLen/2, z+height];
V = [V; x+XLen, y-YLen/2, z+height];

F = [4, 3, 1, 2];
F = [F; 8, 7, 5, 6];
F = [F; 6, 5, 1, 2];
F = [F; 8, 7, 3, 4];
F = [F; 8, 6, 2, 4];
F = [F; 7, 5, 1, 3];

figure(2)
subplot(1,2,1);
hold off;
patch('vertices', V, 'faces', F, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .9 );
camproj( 'perspective');
axis equal
axis off
campos([10,10,10])
hold on;


% might be valuable to label facets as boundary or not, perhaps in more
% compliated cases?
filename = '3d/data/divingBoardM3';
fname = '3d/data/divingBoardM3.poly';
%writePOLY_tetgen(fname,V,F,[]);%,'RegionList',R);

% process the poly file!
%system( pwd + "/lib/svgmesh/tetgen.exe -A -p -q30 -a0.006 " + filename + ".poly" );


[TV,I] = readNODE( filename + ".1.node" );
[TE,A] = readELE( filename + ".1.ele");

subplot(1,2,2)
if numel(A) > 0
    colours = [ 0.5 0.5 0.5; 1 1 0.5; 0.5 1 1; 1 0.5 1; 1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5 ];
    faceColor = colours(mod(A,7)+1,:);
    %%%tetramesh(TE,TV,'faceColor',faceColor,'FaceAlpha',0.1);
    %tetramesh(TE,TV,A,'FaceAlpha',0.1);
    scatter3(TV(:,1),TV(:,2),TV(:,3),'.');
    camproj( 'perspective');
    axis equal
    axis off
    campos([10,10,10])


    %patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', 'flat', 'FaceVertexCData', faceColor, 'FaceAlpha', .5, 'EdgeAlpha', .9 );
else
    %patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .9 );

end
axis equal
axis off


