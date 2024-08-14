function [B, dphidx, Dm, DmInv] = computeB3D(V,T)
    % B = computeB(V,T)
    % computes the kinematic mapping B.
    % V #vertices X 3 entries of x,y,z coordinates of vertices.
    % T #elements X 4 entries of vertex numbers tied to elements.
    numt = size(T,1);
    dphidx = zeros(4,3,numt);
    Dm = zeros(3,3,numt);
    DmInv = zeros(3,3,numt);
    ii = zeros(36*size(T,1),1);
    jj = zeros(36*size(T,1),1);
    vals = zeros(36*size(T,1),1);
    for i = 1:numt
        %Fetching the position of eatch vertex
        t = T(i,:);
        p0 = V(t(1),:);
        p1 = V(t(2),:);
        p2 = V(t(3),:);
        p3 = V(t(4),:);
        
        %femdefo's Dm matrix, this is essentially a material frame of the
        %undeformed mesh
        Dm(:,:,i) = [p1-p0;p2-p0;p3-p0]';
        DmInv(:,:,i) = inv(Dm(:,:,i));
        
        %local B block for individual elements
        D = [-sum(DmInv(:,:,i),1);DmInv(:,:,i)] ;
        dphidx(:,:,i) = D;
        
        %plug the local B into a global B matrix that cont
        elementEnd = i*9;
        cols = reshape([t*3-2;t*3-1;t*3],1,[]);
        colsInds = repmat(cols,3,1);
        jj(36*i-35:36*i) = colsInds;
        
        elementRowInds = repmat([elementEnd-8:elementEnd]',1,4);
        ii(36*i-35:36*i) = elementRowInds(:);
        
        Drep = repmat(D,1,3);
        Dvals = Drep';
        vals(36*i-35:36*i) = Dvals(:);

%         localB = zeros(9,12);
%         localB(1:3,1:3:end) = D';
%         localB(4:6,2:3:end) = D';
%         localB(7:9,3:3:end) = D';
%         B(elementEnd - 8: elementEnd,cols) = localB;
    end
    B = sparse(ii,jj,vals,9*numt, numel(V));
end

