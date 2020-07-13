iter = 150;
close all
for i = 1:iter
    TOPos = [rand*100-50 rand*100-50 rand*20];        % T/O GPS position (including T/O field elevation)
    currPos_rel = [0 0 2]; % Relative GPS NED to T/O position (after auto T/O)
    currPos = [];           % GPS position of drone
    wtPos = [0 0];      % Absolute GPS position of WT
    wtHeight = 70;          % WT nose height from 'ground'
    wtElev = 10;            % WT installation field elevation
    wtAngle = 70;           % yaw angle in degrees
    bladeLen = 30;          % WT blade length
    noFlyWidth = 5;
    WT_path_finding
end