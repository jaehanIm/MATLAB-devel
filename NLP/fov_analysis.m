clear all
% close all
d2r = pi/180;
r2d = 1/d2r;
inspectionDistance = 12;
width = 6;
height = 6 * 1.3143/1.7542;

spdpath = '/home/jaehan/Desktop/210622FT/210622_170128_Sortie7/root_210622_170128.txt';
% gdpath = '/home/jaehan/Desktop/210622FT/210622_170128_Sortie7/gdLog_210622_170128.csv';

% ftdata = readtable(gdpath);

fileID = fopen(spdpath);
% [ftdata,time] = loader(gdpath);

n = 1;
while ~feof(fileID)
    data{n} = fgetl(fileID);
    n = n+1;
end

%% Extraction
posNed = [];
gimbalAngle = [];
yawAngle = [];
shootCmdDuration = [];
compensationAngle = [];

for i = 1:size(data,2)
    % find where it is
    is_posNed = strfind(data{i},'[JobWtMotionBf] Pforplot');
    is_gimbal = strfind(data{i},'[JobWtMotionBf] Gforplot');
    is_yawang = strfind(data{i},'[JobWtMotionBf] Yforplot');
    is_shootCmdDur = strfind(data{i},'shootCmdDuration');
    is_compensation = strfind(data{i},'gimbal compensation angle calculated!');

    
    if ~isempty(is_posNed)
        posNed = vertcat(posNed,str2num(data{i}(is_posNed+25:end)));
    end
    if ~isempty(is_gimbal)
        gimbalAngle = vertcat(gimbalAngle,str2num(data{i}(is_gimbal+25:end)));
    end
    if ~isempty(is_yawang)
        yawAngle = vertcat(yawAngle,str2num(data{i}(is_yawang+25:end)));
    end
    if ~isempty(is_shootCmdDur)
        shootCmdDuration = vertcat(shootCmdDuration,str2num(data{i}(is_shootCmdDur+19:end)));
    end
    if ~isempty(is_compensation)
        compensationAngle = vertcat(compensationAngle,str2num(data{i}(is_compensation+75:end)));
    end
end

compensationAngle = compensationAngle(:,3);

%% Post processing
% autoStartIndex = find(ftdata.fcMcMode == 2 & abs([0;diff(ftdata.fcMcMode)]) == 1);
% compensationAngleTmTc = zeros(size(ftdata,1),1);
% for i = 1:size(autoStartIndex,1)-1
%     compensationAngleTmTc(autoStartIndex(i):autoStartIndex(i+1)-1) = compensationAngle(i);
% end
% compensationAngleTmTc(autoStartIndex(end):end) = compensationAngle(end);

dataLen = size(posNed,1);
zerosData = zeros(dataLen,1);
yawAngleDeg = yawAngle * r2d;
gimbalAngleDeg = gimbalAngle * r2d;

rectangle = [0 0 0 0 0 ; ...
    width/2 width/2 -width/2 -width/2 width/2 ; ...
    height/2 -height/2 -height/2 height/2 height/2];

rectangle(1,:) = inspectionDistance;

dcmI2Hdg = angle2dcm(yawAngle,zerosData,zerosData,'zyx');
dcmHdg2Gimbal = angle2dcm(gimbalAngle(:,3),gimbalAngle(:,1),gimbalAngle(:,2),'zxy');
% dcmHdg2Gimbal = angle2dcm(gimbalAngle(:,3),gimbalAngle(:,2),gimbalAngle(:,1),'zyx'); % How it is shown in GAZEBO

for i = 1:dataLen
    rectHdg(:,:,i) = dcmI2Hdg(:,:,i)' * rectangle;
%     rectNed(:,:,i) = dcmHdg2Gimbal(:,:,i)' * rectHdg(:,:,i);
    rectNed(:,:,i) = dcmI2Hdg(:,:,i)' * dcmHdg2Gimbal(:,:,i)' * rectangle;
end

for i = 1:dataLen
    rectNed(:,:,i) = rectNed(:,:,i)+ posNed(i,:)';
end

%% Plotting
figure(1)
clf
grid on
hold on
plot3(posNed(:,2),posNed(:,1),-posNed(:,3),'r--*')
for i = 1:dataLen
    plot3(rectNed(2,:,i),rectNed(1,:,i),-rectNed(3,:,i),'k')
end
% plot3(ftdata.posNED_m_1,ftdata.posNED_m_0,-ftdata.posNED_m_2,'b');
% xlim([-100 -20])
% ylim([100 160])
% zlim([60 110])

axis equal
view(45,45)

figure(2)
clf
grid on
plot(shootCmdDuration)

figure(3)
clf
grid on
hold on
plot(gimbalAngleDeg(:,1))
plot(gimbalAngleDeg(:,2))
plot(gimbalAngleDeg(:,3))
legend('Roll','Pitch','Yaw')

% figure(4)
% clf
% grid on
% hold on
% plot(ftdata.gimbalRpy_deg_0,'r')
% plot(ftdata.gimbalRpyCmd_deg_0,'r--')
% plot(ftdata.gimbalRpy_deg_1,'b')
% plot(ftdata.gimbalRpyCmd_deg_1,'b--')
% plot(ftdata.gimbalRpy_deg_2 - ftdata.rpy_deg_2 - compensationAngleTmTc,'m')
% plot(ftdata.gimbalRpyCmd_deg_2,'k--')
% ylim([-80 50])