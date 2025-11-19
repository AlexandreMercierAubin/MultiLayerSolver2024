% Make a tube with circular cross section

nSides = 20;
r1 = 1;
r2 = 1.3;
L = 10;

dtheta = 2*pi / nSides;
theta = 0:dtheta:dtheta*(nSides-1);

s = sin(theta);
c = cos(theta);

V =             [ r1*c', r1*s', -L*ones(nSides,1) ];
V = vertcat( V, [ r2*c', r2*s', -L*ones(nSides,1) ] );
V = vertcat( V, [ r1*c', r1*s',  L*ones(nSides,1) ] );
V = vertcat( V, [ r2*c', r2*s',  L*ones(nSides,1) ] );

F = [];
facets = [ 1, 2, nSides+2, nSides+1 ];
for i = 1:nSides-1
    F = vertcat( F, facets+i-1 );
end
F = vertcat( F, [ nSides, 1, nSides+1, nSides*2 ] );

top = F + nSides*2;

F = vertcat( F, top );

% do the sides
side = [ 1, 2, 2 + nSides*2, 1+ nSides*2];
for i = 1:nSides-1
    F = vertcat( F, side  + i - 1 );
    F = vertcat( F, side + nSides + i - 1 );
end
F = vertcat( F, [ nSides, 1, 1 + nSides*2, nSides + nSides*2] );
F = vertcat( F, [ nSides, 1, 1 + nSides*2, nSides + nSides*2] + nSides );



% % Hmm... easier strategy, make the top and bottom and two long strips, then
% % connect top and bottom with various segment strips.  Finally specify the
% % regions.  Triangle doesn't assum CCW faces, so should be fine.
% 
% nR1 = 3;
% nR2 = nR1 - 1;
% width = 1; % width per block
% height = 1;
% totalWidth = (nR1 + nR2) * width;
% n = floor( totalWidth ) + 1;
% 
% z = linspace( 0, totalWidth, n )';
% x0 = zeros( n, 1 );
% x1 = ones( n, 1 ) * height;
% y0 = zeros( n, 1 );
% y1 = ones( n, 1 ) * height;
% V = [ x0 y0 z ];
% V = cat( 1, V, [x1 y0 z] );
% V = cat( 1, V, [x1 y1 z] );
% V = cat( 1, V, [x0 y1 z] );
% 
% facets0 = [ (1:n-1)' (2:n)' (n+2:2*n)' (n+1:2*n-1)' ];
% facets1 = facets0 + n;
% facets2 = facets1 + n;
% facets3 = [ (3*n+1:4*n-1)' (3*n+2:4*n)' (2:n)' (1:n-1)' ];
% 
% F = facets0;
% F = cat( 1, F, facets1 );
% F = cat( 1, F, facets2 );
% F = cat( 1, F, facets3 );
% 
% % do the caps on each end of the cantilever
% 
% cap0 = [ 1, n+1, 2*n+1, 3*n+1 ];
% cap1 = cap0 + n - 1;
% caps = cap0
% caps = cat( 1, caps, cap1 );
% 
% % then add internal boundaries
% 
% for i = 1:n-2
%     cap = cap0 + i;
%     caps = cat( 1, caps, cap );
% end
% 
% F = cat( 1, F, caps );
% 
% % regions
% R = zeros( n-1, 3 );
% R( :, 1 ) = ones(n-1,1)*height/2;
% R( :, 2 ) = ones(n-1,1)*height/2;
% R( :, 3 ) = (1:n-1)'*width - width/2;
% 
% %R( :, 4 ) = (1:n-1)';
% R( 1:2:end, 4 ) = 1;
% R( 2:2:end, 4 ) = 2;


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
filename = '3d/data/tube3dHD';
fname = '3d/data/tube3dHD.poly';
writePOLY_tetgen(fname,V,F,[]);%,'RegionList',R);

% process the poly file!
system( pwd + "/lib/svgmesh/tetgen.exe -A -p -q30 -a0.006 " + filename + ".poly" );


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


