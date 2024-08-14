function rotatedp = rotatepoint(quat,rList)
    %https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html
    q = quaternion(quat);
    qinv = quaternion([quat(1),-quat(2:4)]);
    rQuat = quaternion([zeros(size(rList,1),1),rList]);
    %sometimes people seem to use q^T, but this seems more standard.
    %Changing the order for qpq^-1 yields an active rotation instead of an
    %passive rotation
    rotatedpQuatForm = compact(q*rQuat*qinv);
    rotatedp = rotatedpQuatForm(:,2:4);
end

