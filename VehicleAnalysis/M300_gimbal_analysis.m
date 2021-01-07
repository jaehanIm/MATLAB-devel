gdLog = "/home/jaehan/log/gdLogCsv/gdLog_201022_172206.csv"; % step response
% gdLog = "/home/jaehan/log/gdLogCsv/gdLog_201026_160333.csv"; % smoothGimbal
% gdLog = "/home/jaehan/log/gdLogCsv/gdLog_201028_153059.csv"; % slow smoothGimbal
[data,time]=loader(gdLog);
figure(1)
clf
plot(time,data.gimbalRPY_0,'.--')
hold on
plot(time,data.gimbalRPY_1,'.--')
plot(time,data.gimbalRPY_2,'.--')
grid on
title('Insane Gimbal','fontsize',14)
legend('Roll','Pitch','Yaw')

%% Logging
yawStartFlag = [3467 3987]-9;
yawEndFlag = yawStartFlag + 1/0.02;
pitchStartFlag = [2394 2936]-9;
pitchEndFlag = pitchStartFlag + 1/0.02;
coupleStartFlag = [4500 4946]-9;
coupleEndFlag = coupleStartFlag + 1/0.02;

yawResponse = data.gimbalRPY_2(yawStartFlag(1):yawEndFlag(2));
yawCmd = [zeros(10,1);ones(yawStartFlag(2)-yawStartFlag(1),1);zeros(yawEndFlag(2)-yawStartFlag(2)-9,1)]*30;
pitchResponse = data.gimbalRPY_1(pitchStartFlag(1):pitchEndFlag(2));
pitchCmd = [zeros(10,1);ones(pitchStartFlag(2)-pitchStartFlag(1),1);zeros(pitchEndFlag(2)-pitchStartFlag(2)-9,1)]*(-30);
coupleYawResponse = data.gimbalRPY_1(coupleStartFlag(1):coupleEndFlag(2));
coupleYawCmd = [zeros(10,1);ones(coupleStartFlag(2)-coupleStartFlag(1),1);zeros(coupleEndFlag(2)-coupleStartFlag(2)-9,1)]*(-30);
couplePitchResponse = data.gimbalRPY_2(coupleStartFlag(1):coupleEndFlag(2));
couplePitchCmd = [zeros(10,1);ones(coupleStartFlag(2)-coupleStartFlag(1),1);zeros(coupleEndFlag(2)-coupleStartFlag(2)-9,1)]*30;

rollResponseYaw = data.gimbalRPY_0(yawStartFlag(1):yawEndFlag(2));
rollResponsePitch = data.gimbalRPY_0(pitchStartFlag(1):pitchEndFlag(2));
pitchResponseYaw = data.gimbalRPY_1(yawStartFlag(1):yawEndFlag(2));
yawResponsePitch = data.gimbalRPY_2(pitchStartFlag(1):pitchEndFlag(2));

% detrending
% yawResponse = detrend(yawResponse,0);
% yawCmd = detrend(yawCmd,0);
slip = yawCmd(1); yawCmd = yawCmd - slip;
yawResponse = yawResponse - slip;
slip = pitchCmd(1); pitchCmd = pitchCmd - slip;
pitchResponse = pitchResponse - slip;
slip = coupleYawCmd(1); coupleYawCmd = coupleYawCmd - slip;
coupleYawResponse = coupleYawResponse - slip;
slip = couplePitchCmd(1); couplePitchCmd = couplePitchCmd - slip;
couplePitchResponse = couplePitchResponse - slip;

% tfestimation
tfYaw = tfest(iddata(yawResponse, yawCmd, 0.02),2,0,nan);
tfPitch = tfest(iddata(pitchResponse, pitchCmd, 0.02),2,0,nan);
tfCoupleYaw = tfest(iddata(coupleYawResponse, coupleYawCmd, 0.02),2,0,nan);
tfCouplePitch = tfest(iddata(couplePitchResponse, couplePitchCmd, 0.02),2,0,nan);

% cross tf estimation
tfCrossYawPitch = tfest(iddata(pitchResponseYaw,yawCmd,0.02),2,0,nan);
tfCrossPitchYaw = tfest(iddata(yawResponsePitch, pitchCmd,0.02),2,0,nan);
tfCrossYawRoll = tfest(iddata(rollResponseYaw,yawCmd,0.02),2,0,nan);
tfCrossPitchRoll = tfest(iddata(rollResponsePitch,pitchCmd,0.02),2,0,nan);

%% Plotting

% step response
figure(2)
clf
hold on
grid on
step(tfPitch,'b')
step(tfYaw,'r')
step(tfCouplePitch,'b--')
step(tfCoupleYaw,'r--')
legend('Pitch','Yaw','CouplePitch','CoupleYaw')

figure(3)
clf
hold on
grid on
bode(tfPitch,'b')
bode(tfYaw,'r')
bode(tfCouplePitch,'b--')
bode(tfCoupleYaw,'r--')
legend('Pitch','Yaw','CouplePitch','CoupleYaw')

figure(4)
clf
hold on
grid on
step(tfCrossPitchYaw)
step(tfCrossYawPitch)
step(tfCrossPitchRoll)
step(tfCrossYawRoll)
legend('P->Y','Y->P','P->R','Y->R')

figure(5)
clf
hold on
grid on
rlocus(tfPitch,'b')
rlocus(tfYaw,'r')
rlocus(tfCouplePitch,'b--')
rlocus(tfCoupleYaw,'r--')
legend('Pitch','Yaw','CouplePitch','CoupleYaw')


%% spdLog

% path = convertStringsToChars(path);
% fileID = fopen(path);
% 
% n = 1;
% while ~feof(fileID)
%     data{n} = fgetl(fileID);
%     n = n+1;
% end
% 
% for i = 1:size(data,2)
%     temp = strfind(data{i},'Within');
%     if isempty(temp)
%         remem(i) = nan;
%     else
%         remem(i) = temp;
%     end
% end
% 
% n=1;
% for i = 1:size(data,2)
%     if ~isnan(remem(i))
%         tipDist(n) = str2num(data{i}(76:end-1));
%         n = n+1;
%     end
% end
