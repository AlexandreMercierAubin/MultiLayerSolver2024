function [tanHalfTheta] = tanHalfAngle(n1,n2,n3)
%TANHALFANGLE
    normSub = norm(n1-n2);
    normAdd = norm(n1+n2);
    
    if normAdd == 0
        normComponent = 0;
    else
        normComponent = normSub/normAdd;
    end
    tanHalfTheta = sign(det([n1;n2;n3]'))*normComponent;
end

