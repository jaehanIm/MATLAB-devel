

% Computes Rotation matrix from Euler angle set.
% six basic sequence of rotations around X(Roll),Y(Pitch) and Z(Yaw) axis. 
% Allowed rotations sequences: 
% xyz, xzy, yxz, yzx, zxy, zyx OR rpy, ryp, pry, pyr, yrp, ypr

function DCM=eul2rotm(eul,order)

roll = eul(:,1); pitch = eul(:,2); yaw = eul(:,3);
l = length(roll);
DCM = zeros(3,3,l);

for i = 1:l
    rM=[1    0           0;
        0    cos(roll(i))   sin(roll(i));
        0    -sin(roll(i))  cos(roll(i))];

    pM=[cos(pitch(i))      0       -sin(pitch(i));
        0               1       0;
        sin(pitch(i))      0       cos(pitch(i))];

    yM=[cos(yaw(i))        sin(yaw(i))        0;
        -sin(yaw(i))       cos(yaw(i))        0;
        0               0               1];

    if(strcmpi(order,'xyz')==1 || strcmpi(order,'rpy')==1)
        DCM(:,:,i)=(yM*pM)*rM;  
    elseif(strcmpi(order,'xzy')==1 || strcmpi(order,'ryp')==1)
        DCM(:,:,i)=(pM*yM)*rM;
    elseif(strcmpi(order,'yxz')==1 || strcmpi(order,'pry')==1)
        DCM(:,:,i)=(yM*rM)*pM;
    elseif(strcmpi(order,'yzx')==1 || strcmpi(order,'pyr')==1)
        DCM(:,:,i)=(rM*yM)*pM;
    elseif(strcmpi(order,'zxy')==1 || strcmpi(order,'yrp')==1)
        DCM(:,:,i)=(pM*rM)*yM;
    elseif(strcmpi(order,'zyx')==1 || strcmpi(order,'ypr')==1)
        DCM(:,:,i)=(rM*pM)*yM;
    else
        error('could not recognized the sequence of rotation')
    end
    DCM(:,:,i) = DCM(:,:,i)'; % slight difference to angle2dcm
end
end

