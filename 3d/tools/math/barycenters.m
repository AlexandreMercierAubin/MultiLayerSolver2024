function BC = barycenters( V, T )
    %BC = barycenters( V, T )
    %V: list of 3X#vertices 
    %T: list of triangle vertex indices
    BC = zeros( size(T) );    
    for i = 1:size(T,1)
        BC(i,:) = sum( V(T(i,:),:,1));
    end
end
