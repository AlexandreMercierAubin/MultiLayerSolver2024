function c = cross2D(v1,v2)
%2D cross product, also known as the perp dot product, exterior product or outer product.
    assert(size(v1,2)==2);
    assert(size(v2,2)==2);
    leftRightProducts = v1.*v2(:,[2,1],:);
    c = leftRightProducts(:,1,:)-leftRightProducts(:,2,:);
end
