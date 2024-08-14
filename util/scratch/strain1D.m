function [springStrain] = strain1D(x,springs,restLength)
    %x 3byN 
    % restLength 3 by 1
    springLength = vecnorm(x(springs(:,1)) - x(springs(:,2)),2,1);
    springStrain = (springLength-restLength)./restLength;
end

