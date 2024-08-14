function mesh = generateCantilever2D(attributes, materials, precision, scale, distMeshDimensions, meshingDirection)
    %generateCloth(attributes, materials, precision, scale)
    %generates some 2D rectangular mesh and project it to 3D 
    if nargin < 4
        scale = [1,1,1];
    end
    if nargin < 5
        distMeshDimensions = [3,1]; %size in x and y
    end
    if nargin <6
        meshingDirection = 1;
    end
    halfdim = distMeshDimensions./2;
    %flip halfDim to have the edge line up
    if meshingDirection ==2
        tmp = halfdim(1);
        halfdim(1) = halfdim(2);
        halfdim(2) = tmp;
    end

%     fd = @(p) drectangle(p, -halfdim(1), halfdim(1), -halfdim(2), halfdim(2));
     fd = @(p) drectangle0(p, -halfdim(1), halfdim(1), -halfdim(2), halfdim(2));
    [p, T] = distmesh2d(fd, @huniform, precision, [-halfdim(1), -halfdim(2); halfdim(1), halfdim(2)], [-halfdim(1), -halfdim(2); halfdim(1), -halfdim(2); -halfdim(1), halfdim(2); halfdim(1), halfdim(2)]);
    
    if meshingDirection ==2
        tmp = p(:,1);
        p(:,1) = p(:,2);
        p(:,2)= tmp;
    end

    linePrecision = 1e-6;
    edge1 = p(T(:,2),:) - p(T(:,1),:);
    edge2 = p(T(:,3),:) - p(T(:,1),:);
    edge3 = p(T(:,3),:) - p(T(:,2),:);
    sizeEdge1 = normrow(edge1);
    sizeEdge2 = normrow(edge2);
    sizeEdge3 = normrow(edge3);
    isLine1 = sizeEdge1 <= linePrecision ;
    isLine2 = sizeEdge2 <= linePrecision ;
    isLine3 = sizeEdge3 <= linePrecision ;
    badRatio1 = min(sizeEdge1./sizeEdge2,sizeEdge2./sizeEdge1) < 0.1;
    badRatio2 = min(sizeEdge1./sizeEdge3,sizeEdge3./sizeEdge1) < 0.1;
    badRatio3 = min(sizeEdge2./sizeEdge3,sizeEdge3./sizeEdge2) < 0.1;

    %removes elements that are in the shape of a line
    T(isLine1|isLine2|isLine3|badRatio1|badRatio2|badRatio3, :) = [];

    V = [p(:,1),p(:,2)];

    V(:,1) = scale(1)*V(:,1);
    V(:,2) = scale(2)*V(:,2);

    testAngles = internalangles(V,T);
    if ~isreal(testAngles)
        error("the mesh generated has bad angles... try a different distmesh precision");
    end

    J = [1:size(T,1)]';
    mesh = Mesh(V,T, attributes, materials);
end