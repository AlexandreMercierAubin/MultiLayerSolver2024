function debugDot(p,color,size)
    if nargin < 3
        size = 20;
    end
    plot(p(1:2:end),p(2:2:end),'r.','MarkerSize',size, 'Color',color);
end