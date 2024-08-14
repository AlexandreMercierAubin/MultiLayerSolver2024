function [B, dphidx, Dm] = computeBTri3D(V,T)
    %COMPUTEB Slow but only precomputation. Could probably be vectorized.
    %Not a huge bottleneck for gigantic meshes (adjagency matrix is worse)
   numt = size(T,1);
   dphidx = zeros(4,3,numt);
   Dm = zeros(2,2,numt);

   ii = zeros(27*size(T,1),1);
   jj = zeros(27*size(T,1),1);
   vals = zeros(27*size(T,1),1);
   for i = 1:numt
        %Fetching the position of eatch vertex
        t = T(i,:);
        p0 = V(t(1),:);
        p1 = V(t(2),:);
        p2 = V(t(3),:);
        
        %femdefo's Dm matrix, this is essentially a material frame of the
        %undeformed mesh
        localT = [p1-p0;p2-p0]';
        Dm(:,:,i) = localT'*localT;
        DmInvLocalT = Dm(:,:,i)\localT';
        
        %local B block for individual elements
        D = [-sum(DmInvLocalT,1);DmInvLocalT] ;
        dphidx(:,:,i) = [D;zeros(1,3)]; %adds a row of zeros so it can be stacked with tet D 

        %plug the local B into a global B matrix that cont

        elementEnd = i*9;
        cols = reshape([t*3-2;t*3-1;t*3],1,[]);
        colsInds = repmat(cols,3,1);
        jj(27*i-26:27*i) = colsInds;
        
        elementRowInds = repmat([elementEnd-8:elementEnd]',1,3);
        ii(27*i-26:27*i) = elementRowInds(:);
        
        Drep = repmat(D,1,3);
        Dvals = Drep';
        vals(27*i-26:27*i) = Dvals(:);

%         localB = zeros(9,9);
%         localB(1:3,1:3:end) = D';
%         localB(4:6,2:3:end) = D';
%         localB(7:9,3:3:end) = D';
        
%         rowEnd = i*9;
%         cols = reshape([t*3-2;t*3-1;t*3],1,[]);
%         B(rowEnd - 8: rowEnd,cols) = localB;
   end
   B = sparse(ii,jj,vals,9*numt, numel(V));
end