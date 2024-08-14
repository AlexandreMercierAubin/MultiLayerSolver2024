function debugPlotvectors(p,v,color, useQuivers)
    P = reshape(p, 2, [])';
    V = reshape(v, 2, [])';
    if nargin >= 4 && useQuivers 
        quiver(P(:,1),P(:,2),V(:,1),V(:,2),"Color",color);
    else
        for plotID = 1:size(P,1)
            plot([P(plotID,1),P(plotID,1)+V(plotID,1)],[P(plotID,2),P(plotID,2)+V(plotID,2)],"Color",color);
        end
    end
end