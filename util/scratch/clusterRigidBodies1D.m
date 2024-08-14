function [vertexRigidBodyID, rigidSprings, counter] = clusterRigidBodies1D(strainRates,springs,numRigidSprings)
    [sortedStrains, indices]=sort(strainRates);
    strainRateLimit = sortedStrains(numRigidSprings);
    rigidSprings = find(strainRates <= strainRateLimit);

    vertexRigidBodyID = zeros(max(springs(:)),1);
    counter = 0;
    
    for i = 1 : numel(rigidSprings)%gross quadratic algo for quick prototyping
        springID = rigidSprings(i);
        vert1 = springs(springID,1);
        vert2 = springs(springID,2);

        if(vertexRigidBodyID(vert1) == 0 && vertexRigidBodyID(vert2) == 0)
            counter = counter + 1;
            vertexRigidBodyID(vert1) = counter;
            vertexRigidBodyID(vert2) = counter;
        elseif vertexRigidBodyID(vert1) ~= 0
            vertexRigidBodyID(vert2) = vertexRigidBodyID(vert1);
        elseif vertexRigidBodyID(vert2) ~= 0
            vertexRigidBodyID(vert1) = vertexRigidBodyID(vert2);
        end

        for j = i : numel(rigidSprings)
            springID = rigidSprings(j);
            vert1 = springs(springID,1);
            vert2 = springs(springID,2);
    
            if vertexRigidBodyID(vert1) ~= 0
                vertexRigidBodyID(vert2) = vertexRigidBodyID(vert1);
            elseif vertexRigidBodyID(vert2) ~= 0
                vertexRigidBodyID(vert1) = vertexRigidBodyID(vert2);
            end
        end
    end
end

