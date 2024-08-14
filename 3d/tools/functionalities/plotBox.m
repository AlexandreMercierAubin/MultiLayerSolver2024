function [handle] = plotBox(center, radius, isTriangular, handle)
    p1 = center - radius;
    p2 = p1;
    p2(2) = p2(2) + 2*radius(2);
    p3 = p1;
    p3(1) = p3(1) + 2*radius(1);
    p4 = p1;
    p4(1:2) = p4(1:2) + 2*radius(1:2);
    p5 = p1;
    p5(3) = p5(3) + 2*radius(3);
    p6 = p1;
    p6([2,3]) = p6([2,3]) + 2*radius([2,3]);
    p7 = p1;
    p7([1,3]) = p7([1,3]) + 2*radius([1,3]);
    p8 = center + radius;
    nV = [p1;p2;p3;p4;p5;p6;p7;p8];

    %use this if you want to plot triangles instead;
    if isTriangular
        faces = [1,3,4;
                 1,4,2;
                 1,2,6;
                 1,6,5;
                 1,5,7;
                 1,7,3;
                 2,4,6;
                 4,8,6;
                 3,4,8;
                 8,3,7;
                 5,6,7;
                 6,8,7];
    else
        faces = [1,2,6,5;
                 5,6,8,7;
                 1,2,4,3;
                 1,5,7,3;
                 2,4,8,6;
                 4,8,7,3];
    end

    if nargin < 4 | ~isgraphics(handle)
        hold on;
        handle = patch('vertices', nV, 'faces', faces, 'FaceAlpha', 1, 'EdgeColor', 'red','FaceColor','none','LineWidth',1);
    else
        handle.Vertices = nV;
    end
end

