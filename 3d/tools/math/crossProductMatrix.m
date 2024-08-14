function [crossMat] = crossProductMatrix(m)
% makes a cross product matrix from a vector omega
assert(size(m,2)==3);
z = zeros(1,1,size(m,3));
crossMat = [z,   -m(1,3,:),   m(1,2,:);
            m(1,3,:),    z,  -m(1,1,:);
           -m(1,2,:) ,m(1,1,:),     z];
end

