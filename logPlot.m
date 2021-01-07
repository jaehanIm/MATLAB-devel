%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GD Log Data Analysis
%
% Lee, Byung-Yoon, NearthLab, 181022
% Last modified date: 200803
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all; 
clc;

d2r = pi/180;
r2d = 180/pi;

% Initialization
% File parameters
gdLogFile = '/home/jaehan/log/201105_111716/gdLogCsv/gdLog_201105_111716.csv'; 
% gdLogFile = '/home/jaehan/Desktop/test flight/gdLog_200702_133252.csv'; % 5-1
% gdLogFile = '/home/jaehan/Desktop/test flight/gdLog_200702_123742.csv'; % 4-1
gdLogImageFolder = '/home/jaehan/log/201105_111716/images';
gdLogFileName = gdLogFile;
gdLogFIleMat = gdLogFile + ".mat";

% Time parameters
% dt = 0.02;
timeFormat = 1; % Time format (1: seconds, 2: minute, 3: date)

% Select figures to show
timeStepVer_show = false;
attitudeResp_show = true;
velocityRespWithEulerAngle_show = true;
velocityRespWithGpsVel_show = true;
positionResp_show = true;
errLat_show = false;
bdLidar_show = false;
acLidarH_show = false;
acLidarV_show = false;
windSensor_show = false;
plot3dEstimatedBlade_show = false;
plot3d_show = false;
plot3dPosCmd_show = false;
situationAnalysis3d_show = false;
remoteController_show = true;
gimbalAngle_show = false;
bladeTravelDistance_show = false;
missionTime_show = false;
analEachJob_show = false;


%%
% Load Log Data
gdLog = importfile(gdLogFileName);
data = gdLog;
timeValidIdx = max(find(gdLog.rosTime < 1514764800)); % 2018/01/01 == 1514764800
data((1:timeValidIdx),:) = [];
dataSize = size(data.rosTime);

fTimeGlobal = data.rosTime;
fTimeGlobal_Date = datetime(data.rosTime,'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
fTimeLocal = fTimeGlobal(:,1) - fTimeGlobal(1,1);
fTimeLocal_Sec = seconds(fTimeLocal);

time = fTimeLocal_Sec;
timeSec = fTimeLocal_Sec;
hourTime = fTimeGlobal_Date;
dateStart = fTimeGlobal_Date(1);

switch timeFormat % (1: seconds, 2: minute, 3: date)
    case 1
        strTimeXtickFormat = 's';
    case 2
        strTimeXtickFormat = 'mm:ss.SSS';
    case 3
        time = time + dateStart;
        strTimeXtickFormat = 'HH:mm:ss';
    otherwise
end

flightMode = data.flightMode;
ctrlDeviceStatus = data.ctrlDeviceStatus;
fcMcMode = data.fcMcMode;
nSat = data.nSat;
gpsFix = data.gpsFix;
jobSeq = data.jobSeq;
ySpType = data.ySpType;
velNedGps(:,1) = data.velNedGps_0;
velNedGps(:,2) = data.velNedGps_1;
velNedGps(:,3) = data.velNedGps_2;
posNed(:,1) = data.posNed_0;
posNed(:,2) = data.posNed_1;
posNed(:,3) = data.posNed_2;
velNed(:,1) = data.velNed_0;
velNed(:,2) = data.velNed_1;
velNed(:,3) = data.velNed_2;
rpy(:,1) = data.rpy_0;
rpy(:,2) = data.rpy_1;
rpy(:,3) = data.rpy_2;
rpyRad = rpy*d2r;
dcmI2B = angle2dcm(rpyRad(:,3),rpyRad(:,2),rpyRad(:,1),'ZYX');
dcmI2N = angle2dcm(rpyRad(:,3), zeros(size(rpyRad,1),1), zeros(size(rpyRad,1),1),'ZYX');
velNav = zeros(size(velNed,1),3);
for i=1:size(velNed,1)
    velNav(i,:)=dcmI2N(:,:,i)*velNed(i,:)';
end
velNavGps = zeros(size(velNedGps,1),3);
for i=1:size(velNedGps,1)
    velNavGps(i,:)=dcmI2N(:,:,i)*velNedGps(i,:)';
end
posNav = zeros(size(posNed,1),3);
for i = 1:size(posNed,1)
    posNav(i,:) = dcmI2N(:,:,i)*posNed(i,:)';
end
ctrlStruct = data.ctrlStruct;
ctrlSetpointType = data.ctrlSetpointType;
ctrlOutputType = data.ctrlOutputType;
ctrlUser = data.ctrlUser;
ctrlSp(:,1) = data.ctrlSp_0;
ctrlSp(:,2) = data.ctrlSp_1;
ctrlSp(:,3) = data.ctrlSp_2;
ctrlOp(:,1) = data.ctrlOp_0;
ctrlOp(:,2) = data.ctrlOp_1;
ctrlOp(:,3) = data.ctrlOp_2;
ySp = data.ySp;
rcRoll = data.rcRoll;
rcPitch = data.rcPitch;
rcYaw = data.rcYaw;
rcThrottle = data.rcThrottle;
GpsNSV = data.GpsNSV;
RtkHealthFlag = data.RtkHealthFlag;
GpsFusedNSV = data.GpsFusedNSV;
GpHealth = data.GpHealth;
posGPS(:,1) = data.posGPS_0*10^-7;
posGPS(:,2) = data.posGPS_1*10^-7;
posGPS(:,3) = data.posGPS_2*10^-3;
posRTK(:,1) = data.posRTK_0;
posRTK(:,2) = data.posRTK_1;
posRTK(:,3) = data.posRTK_2;
posGpsFused(:,1) = data.posGpsFused_0*r2d;
posGpsFused(:,2) = data.posGpsFused_1*r2d;
posGpsFused(:,3) = data.posGpsFused_2;
posGP(:,1) = data.posGp_0*r2d;
posGP(:,2) = data.posGp_1*r2d;
posGP(:,3) = data.posGp_2;
velUvw(:,1) = data.vbx;
velUvw(:,2) = data.vby;
velUvw(:,3) = data.vbz;
acWarnStat = data.AcWarnStat;
acHorRange = data.AcHorWarnRange;
acHorAngle = data.AcHorWarnAngle;
acVerRange = data.AcVerWarnRange;
acVerAngle = data.AcVerWarnAngle;
bdLidarDist = data.LidarDist;
bdLidarAngle = data.LidarAngle;
errLatMix = data.errLatMix;
errLatVis = data.errLatVis;
errLatLid = data.errLatLid;
cmdLatVelIgain = data.cmdLatVelIgain;
cmdLatVelMix = data.cmdLatVelMix;
velCtrlI_u = data.velCtrlI_u;
velCtrlI_v = data.velCtrlI_v;
posCtrlI(:,1) = data.posCtrlI_N;
posCtrlI(:,2) = data.posCtrlI_E;
posCtrlI(:,3) = data.posCtrlI_D;
gimbalRPY(:,1) = data.gimbalRPY_0;
gimbalRPY(:,2) = data.gimbalRPY_1;
gimbalRPY(:,3) = data.gimbalRPY_2;
windStatus = data.windStatus;
windSpeed = data.windSpeed;
windAngle = data.windAngle;
windQueryTime = data.windQueryTime;
windResponseTime = data.windResponseTime;
acousticTemp = data.acousticTemp;
tempQueryTime = data.tempQueryTime;
tempResponseTime = data.tempResponseTime;
accBody(:,1) = data.accBody_0;
accBody(:,2) = data.accBody_1;
accBody(:,3) = data.accBody_2;
trajUnitVectorT(:,1) = data.trajUnitVectorT_0;
trajUnitVectorT(:,2) = data.trajUnitVectorT_1;
trajUnitVectorT(:,3) = data.trajUnitVectorT_2;
trajUnitVectorN(:,1) = data.trajUnitVectorN_0;
trajUnitVectorN(:,2) = data.trajUnitVectorN_1;
trajUnitVectorN(:,3) = data.trajUnitVectorN_2;
trajUnitVectorB(:,1) = data.trajUnitVectorB_0;
trajUnitVectorB(:,2) = data.trajUnitVectorB_1;
trajUnitVectorB(:,3) = data.trajUnitVectorB_2;
trajCmdTnb(:,1) = data.trajCtrlI_T;
trajCmdTnb(:,2) = data.trajCtrlI_N;
trajCmdTnb(:,3) = data.trajCtrlI_B;
StdJobLongPidErr = data.StdJobLongPidErr;
StdJobLongPidRate = data.StdJobLongPidRate;
StdJobLongPidIgain = data.StdJobLongPidIgain;
GuideModeLongPidErr = data.GuideModeLongPidErr;
GuideModeLongPidRate = data.GuideModeLongPidRate;
GuideModeLongPidIgain = data.GuideModeLongPidIgain;
pqr(:,1) = data.pqr_0;
pqr(:,2) = data.pqr_1;
pqr(:,3) = data.pqr_2;
rpdCmd(:,1) = data.rpdCmd_0;
rpdCmd(:,2) = data.rpdCmd_1;
rpdCmd(:,3) = data.rpdCmd_2;
velCmdNav(:,1) = data.velCmdNav_0;
velCmdNav(:,2) = data.velCmdNav_1;
velCmdNav(:,3) = data.velCmdNav_2;
posCmdNed(:,1) = data.posCmdNed_0;
posCmdNed(:,2) = data.posCmdNed_1;
posCmdNed(:,3) = data.posCmdNed_2;
trajTimeCur = data.trajTimeCur;
trajTimeMax = data.trajTimeMax;
bladeTravelDistance = data.bladeTravelDistance;

sTime = min(time);
eTime = max(time);
timeDurTotal = [sTime eTime];

%%
% Time Step Verification
if timeStepVer_show
    timeStep = zeros(dataSize-1);
    for i = 1:(dataSize-1)
        timeStep(i) = data.rosTime(i+1) - data.rosTime(i);
    end
    timeStepMax = max(timeStep);
    timeStepMin = min(timeStep);

    figure(1);
    set(gcf, 'Position', [100, -20, 800, 800]);
    subFigTime(1) = subplot(2,1,1);
    hEdges = 0:0.001:0.14;
    histogram(timeStep, hEdges, 'Normalization', 'probability'); 
    ylim([0 0.7]);
    title("FC Time Step");
    xlabel("Time Step (sec)");
    ylabel("Probability");
    timeStepStr = "Desired Operation Rate: 50Hz (0.02 sec)";
    timeStepStr = timeStepStr + newline + "max time step: " + string(timeStepMax) + " (sec)";
    timeStepStr = timeStepStr + newline + "min time step: " + string(timeStepMin) + " (sec)";
    legend(timeStepStr);
    grid on;
    
    subFigTime(2) = subplot(2,1,2);
    plot(time(1:end-1,1), seconds(diff(time))); ylim([0 0.1]);
    xtickformat(strTimeXtickFormat);
    legend("dt");
    grid on;
end

%%
% Attitude Responses
if attitudeResp_show
    figure(2); 
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("Attitude Responses");
    
    subFigAtt(1) = subplot(4,1,1);
    plot(time, jobSeq, '-k', 'LineWidth', 1.5); hold on; 
    plot(time, fcMcMode, '--b', 'LineWidth', 1.5); hold on;
    plot(time, ctrlSetpointType, '-.r', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    ylim([0 22]);
    legend("jobSeq","fcMcMode","ctrlSetpointType");
    ylabel("Index");
    grid on;

    subFigAtt(2) = subplot(4,1,2);
    plot(time, rpdCmd(:,1), ':r', 'LineWidth', 1.5); hold on;
    plot(time, rpy(:,1), '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    ylim([-20 20]);
    legend("Roll Cmd","Roll Angle");
    ylabel("roll (deg)");
    grid on;

    subFigAtt(3) = subplot(4,1,3);
    plot(time, rpdCmd(:,2), ':r', 'LineWidth', 1.5); hold on;
    plot(time, rpy(:,2), '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    ylim([-20 20]);
    legend("Pitch Cmd","Pitch Angle");
    ylabel("pitch (deg)");
    grid on;

    subFigAtt(4) = subplot(4,1,4);
    plot(time, ySp, ':r', 'LineWidth', 1.5); hold on;
    plot(time, rpy(:,3), '-k', 'LineWidth', 1.2); hold on;
    xtickformat(strTimeXtickFormat);
    ylim([-180 180]);
    legend('Yaw Cmd','Yaw Angle');
    ylabel("yaw angle");
    xlabel("time (sec)");
    grid on;

    linkaxes([subFigAtt(1), subFigAtt(2), subFigAtt(3), subFigAtt(4)], 'x');
    xlim([sTime eTime]);
end

%%
% Velocity Command Responses with Euler Angle (Navigation frame)
if velocityRespWithEulerAngle_show
    figure(3);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("Velocity Command Responses with Euler Angle Command");
    
    subFigVelCmd(1) = subplot(3,1,1);
    plot(time, jobSeq, '-k', 'LineWidth', 1.5); hold on; 
    plot(time, fcMcMode, '--b', 'LineWidth', 1.5); hold on;
    plot(time, ctrlSetpointType, '-.r', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    ylim([0 22]);
    legend("jobSeq","fcMcMode","ctrlSetpointType");
    ylabel("Index");
    grid on;

    subFigVelCmd(2) = subplot(3,1,2);
    plot(time, velCmdNav(:,1), ':r', 'LineWidth', 1.5); hold on;
    plot(time, velNav(:,1), '-k', 'LineWidth', 1.2); 
    plot(time, -rpdCmd(:,2)*d2r*3, '-.b', 'LineWidth', 1.2); hold on;
    plot(time, -velCtrlI_u(:,1)*3, '--g', 'LineWidth', 1.5); hold on;
    xtickformat(strTimeXtickFormat);
    ylim([-1.5 1.5]);
    legend("Velocity Sp (NavX)","Velocity (NavX)","-Pitch Cmd(rad)*3","-Pitch Cmd I(Vel.Ctrl. I)*3");
    ylabel("Velocity (Nav X) (m/s)");
    grid on;

    subFigVelCmd(3) = subplot(3,1,3);
    plot(time, velCmdNav(:,2), ':r', 'LineWidth', 1.5); hold on;
    plot(time, velNav(:,2), '-k', 'LineWidth', 1.2); 
    plot(time, rpdCmd(:,1)*d2r*3, '-.b', 'LineWidth', 1.2); hold on;
    plot(time, velCtrlI_v(:,1)*3, '--g', 'LineWidth', 1.5); hold on;
    xtickformat(strTimeXtickFormat);
    ylim([-1.5 1.5]);
    legend("Velocity Sp (NavY)","Velocity (NavY)","Roll Cmd(rad)*3","Roll Cmd I(Vel.Ctrl. I)*3");
    ylabel("Velocity (Nav Y) (m/s)");
    xlabel("time (sec)");
    grid on;

    linkaxes([subFigVelCmd(1), subFigVelCmd(2), subFigVelCmd(3)], 'x');
end

%%
% Velocity Command Responses with GPS Velocity (Navigation frame)
if velocityRespWithGpsVel_show
    figure(4);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("Velocity Command Responses with GPS Velocity");

    subFigVel(1) = subplot(4,1,1);
    plot(time, jobSeq, '-k', 'LineWidth', 1.5); hold on; 
    plot(time, fcMcMode, '--b', 'LineWidth', 1.5); hold on;
    plot(time, ctrlSetpointType, '-.r', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    ylim([0 22]);
    legend("jobSeq","fcMcMode","ctrlSetpointType");
    ylabel("Index");
    grid on;

    subFigVel(2) = subplot(4,1,2);
    plot(time, velCmdNav(:,1), ':r', 'LineWidth', 1.5); hold on;
    plot(time, velNav(:,1), '-k', 'LineWidth', 1.2); hold on;
    plot(time, velNavGps(:,1), '-.b', 'LineWidth', 1.2); 
    xtickformat(strTimeXtickFormat);
    ylim([-1.5 1.5]);
    legend("Velocity Sp (NavX)","Velocity (NavX)","GPS Velocity (NavX)");
    ylabel("Velocity (Nav X) (m/s)");
    grid on;

    subFigVel(3) = subplot(4,1,3);
    plot(time, velCmdNav(:,2), ':r', 'LineWidth', 1.5); hold on;
    plot(time, velNav(:,2), '-k', 'LineWidth', 1.2); hold on;
    plot(time, velNavGps(:,2), '-.b', 'LineWidth', 1.2); 
    xtickformat(strTimeXtickFormat);
    ylim([-1.5 1.5]);
    legend("Velocity Sp (NavY)","Velocity (NavY)","GPS Velocity (NavY)");
    ylabel("Velocity (Nav Y) (m/s)");
    grid on;

    subFigVel(4) = subplot(4,1,4);
    plot(time, velCmdNav(:,3), ':r', 'LineWidth', 1.5); hold on;
    plot(time, velNav(:,3), '-k', 'LineWidth', 1.2); hold on;
    plot(time, velNavGps(:,3), '-.b', 'LineWidth', 1.2); 
    xtickformat(strTimeXtickFormat);
    ylim([-1.5 1.5]);
    legend("Velocity Sp (NavZ)","Velocity (NavZ)","GPS Velocity (NavZ)");
    ylabel("Velocity (Nav Z) (m/s)");
    xlabel("time (sec)");
    grid on;

    linkaxes([subFigVel(1), subFigVel(2), subFigVel(3), subFigVel(4)], 'x');
end

%%
% Position Responses
if positionResp_show
    figure(5);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("Position Responses");
    
    subFigPos(1) = subplot(4,1,1);
    plot(time, jobSeq, '-k', 'LineWidth', 1.5); hold on; 
    plot(time, fcMcMode, '--b', 'LineWidth', 1.5); hold on;
    plot(time, ctrlSetpointType, '-.r', 'LineWidth', 1.5); hold on;
%     plot(time, ctrlUser, '-.g', 'LineWidth', 1.2);
%     plot(time, flightMode, ':m', 'LineWidth', 1.5);
    plot(time, nSat, '-.g', 'LineWidth', 1.2);
    plot(time, gpsFix, ':m', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    ylim([0 22]);
    legend("jobSeq","fcMcMode","ctrlSetpointType", "nSat", "gpsFix");
    ylabel("Index");
    grid on;

    subFigPos(2) = subplot(4,1,2);
    plot(time, posCmdNed(:,1), '--r', 'LineWidth', 1.5); hold on;
    plot(time, posCtrlI(:,1)+posCmdNed(:,1), '--b', 'LineWidth', 1.5); hold on;
    plot(time, posNed(:,1), '-k', 'LineWidth', 1.2); 
    xtickformat(strTimeXtickFormat);
    ylim([-inf inf]);
    legend("Position Sp (N)", "Position Sp+I (N)", "Position (N)");
    ylabel("Position (N) (m)");
    grid on;

    subFigPos(3) = subplot(4,1,3);
    plot(time, posCmdNed(:,2), '--r', 'LineWidth', 1.5); hold on;
    plot(time, posCtrlI(:,2)+posCmdNed(:,2), '--b', 'LineWidth', 1.5); hold on;
    plot(time, posNed(:,2), '-k', 'LineWidth', 1.2); 
    xtickformat(strTimeXtickFormat);
    ylim([-inf inf]);
    legend("Position Sp (E)", "Position Sp+I (E)","Position (E)");
    ylabel("Position (E) (m)");
    grid on;

    subFigPos(4) = subplot(4,1,4);
    plot(time, posCmdNed(:,3), '--r', 'LineWidth', 1.5); hold on;
    plot(time, posCtrlI(:,3)+posCmdNed(:,3), '--b', 'LineWidth', 1.5); hold on;
    plot(time, posNed(:,3), '-k', 'LineWidth', 1.2); 
    xtickformat(strTimeXtickFormat);
    ylim([-inf inf]);
    legend("Position Sp (D)", "Position Sp+I (D)","Position (D)");
    ylabel("Position (D) (m)");
    xlabel("time (sec)");
    grid on;

    linkaxes([subFigPos(1), subFigPos(2), subFigPos(3), subFigPos(4)], 'x');
end

%%
% errLat Responses
if errLat_show
    figure(6);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("errLat Responses");
    
    subFigErrLat(1) = subplot(3,1,1);
    plot(time, jobSeq, '-k', 'LineWidth', 1.5); hold on; 
    plot(time, fcMcMode, '--b', 'LineWidth', 1.5); hold on;
    plot(time, ctrlSetpointType, '-.r', 'LineWidth', 1.5); hold on;
    plot(time, ctrlUser, '--g', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    ylim([0 22]);
    legend("jobSeq","fcMcMode","ctrlSetpointType", "ctrlUser");
    ylabel("Index");
    grid on;

    subFigErrLat(2) = subplot(3,1,2);
    plot(time, errLatVis, '--r', 'LineWidth', 1.5); hold on;
    plot(time, errLatLid, '--b', 'LineWidth', 1.5); hold on;
    plot(time, errLatMix, '-k', 'LineWidth', 1.2); 
    xtickformat(strTimeXtickFormat);
    ylim([-1.5 1.5]);
    legend("errLatVis", "errLatLid", "errLatMix");
    ylabel("errLat");
    grid on;

    subFigErrLat(3) = subplot(3,1,3);
    plot(time, cmdLatVelIgain, '--r', 'LineWidth', 1.5); hold on;
    plot(time, cmdLatVelMix, '--b', 'LineWidth', 1.5); hold on;
    xtickformat(strTimeXtickFormat);
    ylim([-1.5 1.5]);
    legend("cmdLatVelIgain", "cmdLatVelMix");
    ylabel("cmdLatVel");
    grid on;

    linkaxes([subFigErrLat(1), subFigErrLat(2), subFigErrLat(3)], 'x');
end

%%
% BD Lidar
if bdLidar_show
    figure(7);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("BD Lidar Information");

    subFigLid(1) = subplot(2,1,1);
    plot(time, bdLidarDist, '-k', 'LineWidth', 1.2); grid on;
    xtickformat(strTimeXtickFormat);
    legend("BD Lidar Distance");
    ylabel("Distance (m)");
    xlabel("time (sec)");
    
    subFigLid(2) = subplot(2,1,2);
    plot(time, bdLidarAngle, '-k', 'LineWidth', 1.2); grid on;
    xtickformat(strTimeXtickFormat);
    legend("BD Lidar Angle");
    ylabel("Angle (deg)");
    xlabel("time (sec)");

    linkaxes([subFigLid(1), subFigLid(2)], 'x');
end

%%
% AC Hor Lidar
if acLidarH_show
    figure(8);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("AC Hor Lidar Information");

    subFigAcHorLid(1) = subplot(2,1,1);
    plot(time, acHorRange, '-k', 'LineWidth', 1.2); grid on; hold on;
    plot(time, acWarnStat, '--r', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    legend("AC Hor Lidar Distance", "AC Warn Stat");
    ylabel("Distance (m)");
    xlabel("time (sec)");
    
    subFigAcHorLid(2) = subplot(2,1,2);
    plot(time, acHorAngle, '-k', 'LineWidth', 1.2); grid on;
    xtickformat(strTimeXtickFormat);
    legend("AC Hor Lidar Angle");
    ylabel("Angle (deg)");
    xlabel("time (sec)");

    linkaxes([subFigAcHorLid(1), subFigAcHorLid(2)], 'x');
    
%     size(find( (acHorRange(100/dt:550/dt)>2) & (acHorRange(100/dt:550/dt)<4.5)))
end

%%
% AC Ver Lidar
if acLidarV_show
    figure(9);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("AC Ver Lidar Information");

    subFigAcVerLid(1) = subplot(2,1,1);
    plot(time, acVerRange, '-k', 'LineWidth', 1.2); grid on; hold on;
    plot(time, acWarnStat, '--r', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    legend("AC Ver Lidar Distance", "AC Warn Stat");
    ylabel("Distance (m)");
    xlabel("time (sec)");
    
    subFigAcVerLid(2) = subplot(2,1,2);
    plot(time, acVerAngle, '-k', 'LineWidth', 1.2); grid on;
    xtickformat(strTimeXtickFormat);
    legend("AC Ver Lidar Angle");
    ylabel("Angle (deg)");
    xlabel("time (sec)");

    linkaxes([subFigAcVerLid(1), subFigAcVerLid(2)], 'x');
    
%     size(find( (acVerRange(100/dt:550/dt)>2) & (acVerRange(100/dt:550/dt)<4.5)))
end

%% 
% Wind Sensor Related Plots
if windSensor_show
    figure(10);
    
    plot(time, windStatus);
    xtickformat(strTimeXtickFormat);
    xlabel("time (sec)");
    ylabel("Wind Sensor Data Status");
    grid on;
    
    figure(11);
    
    subFigSpeed = subplot(4,1,1);
    
    plot(time, windSpeed);
    xtickformat(strTimeXtickFormat);
    xlabel("time (sec)");
    ylabel("Wind Speed (m/s)");
    grid on;
    
    subFigAngle = subplot(4,1,2);
    
    plot(time, windAngle);
    xtickformat(strTimeXtickFormat);
    xlabel("time (sec)");
    ylabel("Wind Angle (��)");
    grid on;   
    
    subFigTemp = subplot(4,1,3);
    
    plot(time, acousticTemp);
    xtickformat(strTimeXtickFormat);
    xlabel("time (sec)");
    ylabel("Wind Temperature (��C)");
    grid on;

    subFigTime = subplot(4,1,4);
    
    plot(time, windResponseTime - windQueryTime, '-b'); hold on;
    plot(time, tempResponseTime - tempQueryTime, '-r');
    xtickformat(strTimeXtickFormat);
    legend("Wind Time Diff", "Temp Time Diff");
    xlabel("time (sec)");
    ylabel("Wind Query/Response Time Diff(sec)");
    grid on;  
   
    linkaxes([subFigSpeed, subFigAngle, subFigTemp, subFigTime], 'x');
    
    figure(12); % Figure of True Wind Speed and Angle
    % 1. V_Wind = V_wind' + V_Drone + r * bodyrate
    % 2. V_Wind' = [|V_Wind| * cos(windAngle), |V_Wind| * sin(windAngle), 0]
    % 3. V_Drone | proj Sensor Plane = V_Drone - V_Drone | proj plane normal
    
    armLength = 0.3;
    windSensorPlaneNormalInit = [0, 0, -1]';
    windSensorNorthInit = [1, 0, 0]';
    arraySize = size(windSpeed,1);
    
    % 2.
    windVelSensor = [windSpeed .* cos(windAngle .* d2r), windSpeed .* sin(windAngle .* d2r), zeros(arraySize,1)];
    
    % 3.
    rotMatrix = eul2rotm([rpy(:,1), rpy(:,2), rpy(:,3)]*d2r, 'XYZ');
    
    windSensorPlaneNormal = zeros(arraySize, 3);
    windSensorNorth = zeros(arraySize, 3);
    windSensorPlaneNormalNorm = zeros(arraySize, 1);
    dot_NormalDroneVel = zeros(arraySize,1);
    droneBodyVel_ProjSensorPlaneNormal = zeros(arraySize, 3);
    droneBodyVel_ProjSensorPlane = zeros(arraySize, 3);
    windVelFrame = zeros(arraySize, 3);
    
    unsignedWindAngle = zeros(arraySize, 1);
    sgn = zeros(arraySize, 3);
    windAngleTrue = zeros(arraySize, 1);
    
    for i = 1: arraySize
        windSensorNorth(i,:) = transpose(rotMatrix(:,:,i) * windSensorNorthInit);
        windSensorPlaneNormal(i,:) = transpose(rotMatrix(:,:,i) * windSensorPlaneNormalInit);
        dot_NormalDroneVel(i) = dot(windSensorPlaneNormal(i,:), velUvw(i,:));
        windSensorPlaneNormalNorm(i) = sqrt(windSensorPlaneNormal(i,1).^2 + windSensorPlaneNormal(i,2).^2 + windSensorPlaneNormal(i,3).^2);
        droneBodyVel_ProjSensorPlaneNormal(i,:) = dot_NormalDroneVel(i)/windSensorPlaneNormalNorm(i) * windSensorPlaneNormal(i,:);
        droneBodyVel_ProjSensorPlane(i,:) = velUvw(i,:) - droneBodyVel_ProjSensorPlaneNormal(i,:);
        windVelFrame(i,:) = transpose(rotMatrix(:,:,i) * windVelSensor(i,:)');
    end
    
    windVelTrue = windVelFrame - droneBodyVel_ProjSensorPlane + armLength .* pqr;
    windSpeedTrue = sqrt(windVelTrue(:,1).^2 + windVelTrue(:,2).^2 + windVelTrue(:,3).^2);
    minWind = min(windSpeedTrue)
    maxWind = max(windSpeedTrue)
    meanWind = mean(windSpeedTrue)
    % 4. Calculating Angle
    
    for i = 1: arraySize
        unsignedWindAngle(i) = acos(dot(windSensorNorth(i,:), windVelTrue(i,:)/norm(windVelTrue(i,:))));
        sgn(i,:) = sign(cross(windSensorNorth(i,:), windVelTrue(i,:)));
        windAngleTrue(i) = mod(unsignedWindAngle(i) * (-1)^(sgn(i,3)<0), 2*pi);
    end

    subFigTrueSpeedADS = subplot(4,1,1);
    
    plot(time, windSpeedTrue, time, windSpeed, '--r'); hold on;
    ylabel("Wind Speed True(m/s)");
    legend("Wind Speed True", "Wind Speed Sensor");
    grid on;  
    
    subFigTrueAngleADS = subplot(4,1,2);
    
    plot(time, windAngleTrue*r2d, time, windAngle, '--r');hold on;
    xtickformat(strTimeXtickFormat);
    legend("Wind Angle True", "Wind Angle Sensor");
    xlabel("time (sec)");
    ylabel("Wind Angle True(��)");
    grid on;  
    
    subFigVelUVW = subplot(4,1,3);
    
    plot(time, velUvw(:,1), time, velUvw(:,2), time, velUvw(:,3));hold on;
    xtickformat(strTimeXtickFormat);
    legend("U Velocity", "V Velocity", "W Velocity")
    xlabel("time (sec)");
    ylabel("Wind Speed(m/s)");
    grid on;  
    
    subFigPQR = subplot(4,1,4);
    
    plot(time, pqr(:,1)*r2d, time, pqr(:,2)*r2d, time, pqr(:,3)*r2d);hold on;
    xtickformat(strTimeXtickFormat);
    legend("P", "Q", "R")
    xlabel("time (sec)");
    ylabel("PQR(��/s)");
    grid on; 
    
    linkaxes([subFigTrueSpeedADS, subFigTrueAngleADS, subFigVelUVW, subFigPQR], 'x');
end

%% 
% 3d plot with estimated blade
if plot3dEstimatedBlade_show
    figure(13);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("3d Plot with Estimated Blade");

    % for all time
%     sTimeIdx = 1;
%     eTimeIdx = size(time,1);

    % for mission time
    sTimeIdx = min(find(jobSeq==1));
    eTimeIdx = min(find( (jobSeq==max(jobSeq)) & (fcMcMode==0) ));

    % for specific job
%     sTimeIdx = min(find(jobSeq==2));
%     eTimeIdx = max(find(jobSeq==2));
    
    durTimeIdx = [sTimeIdx:eTimeIdx];

    estimatedBladePosNed = zeros(size(posNed,1), 3);

    for i=1:size(bdLidarDist)
        if(bdLidarDist(i,1) > 0)
            bladeRelPosB = [bdLidarDist(i,1); 0; 0];
            lidarRelDcmI2B = angle2dcm((rpy(i,3)*d2r + bdLidarAngle(i,1)*d2r), 0, 0, 'ZYX');
            lidarRelDcmB2I = lidarRelDcmI2B';
            bladeRelPosI = lidarRelDcmB2I*bladeRelPosB;
            estimatedBladePosNed(i,:) = posNed(i,:)' + bladeRelPosI;
        else
            estimatedBladePosNed(i,:) = [0, 0, 0];
        end
    end

    plot3(estimatedBladePosNed([sTimeIdx:eTimeIdx],1), estimatedBladePosNed([sTimeIdx:eTimeIdx],2), -estimatedBladePosNed([sTimeIdx:eTimeIdx],3), 'b.', 'LineWidth', 1.0); grid on; hold on;
    plot3(posNed([sTimeIdx:eTimeIdx],1), posNed([sTimeIdx:eTimeIdx],2), -posNed([sTimeIdx:eTimeIdx],3), '-k', 'LineWidth', 1.5); grid on; hold on;

    % GP
%     posGP_ned = lla2flat(posGP, [posGP(1,1), posGP(2,2)], 0, -posGP(3,3), 'WGS84');
%     plot3(posGP_ned([sTimeIdx:eTimeIdx],1), posGP_ned([sTimeIdx:eTimeIdx],2), -posGP_ned([sTimeIdx:eTimeIdx],3), 'r-'); grid on; hold on;

    % GPS raw
%     posGPS_ned = lla2flat(posGPS, [posGPS(1,1), posGPS(2,2)], 0, -posGPS(3,3), 'WGS84');
%     plot3(posGPS_ned([sTimeIdx:eTimeIdx],1), posGPS_ned([sTimeIdx:eTimeIdx],2), -posGPS_ned([sTimeIdx:eTimeIdx],3), 'g-'); grid on; hold on;

    % RTK raw
%     posRTK_ned = lla2flat(posRTK, [posRTK(1,1), posRTK(2,2)], 0, -posRTK(3,3), 'WGS84');
%     plot3(posRTK_ned([sTimeIdx:eTimeIdx],1), posRTK_ned([sTimeIdx:eTimeIdx],2), -posRTK_ned([sTimeIdx:eTimeIdx],3), 'b-'); grid on; hold on;
    
    xlabel("Position(N)");
    ylabel("Position(E)");
    zlabel("Position(Alt)");
    axis ij;
    axis equal;
    legend('estimated blade','vehicle position');
end


%% 
% 3d plot with position command for trajectory following
if plot3dPosCmd_show
    figure(14);
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("3d Plot with Position Command");

    % for all time
%     sTimeIdx = 1;
%     eTimeIdx = size(time,1);

    % for mission time
%     sTimeIdx = min(find(jobSeq==1));
%     eTimeIdx = min(find( (jobSeq==max(jobSeq)) & (fcMcMode==0) ));

    % for specific job
    sTimeIdx = min(find( (jobSeq==8) & (ctrlSetpointType==10) ));
    eTimeIdx = max(find( (jobSeq==8) & (ctrlSetpointType==10) ));
    durTimeIdx = [sTimeIdx:eTimeIdx];
    
    estimatedBladePosNed = zeros(size(posNed,1), 3);

    for i=1:size(bdLidarDist)
        if(bdLidarDist(i,1) > 0)
            bladeRelPosB = [bdLidarDist(i,1); 0; 0];
            lidarRelDcmI2B = angle2dcm((rpy(i,3)*d2r + bdLidarAngle(i,1)*d2r), 0, 0, 'ZYX');
            lidarRelDcmB2I = lidarRelDcmI2B';
            bladeRelPosI = lidarRelDcmB2I*bladeRelPosB;
            estimatedBladePosNed(i,:) = posNed(i,:)' + bladeRelPosI;
        else
            estimatedBladePosNed(i,:) = [0, 0, 0];
        end
    end

%     plot3(estimatedBladePosNed(durTimeIdx,1), estimatedBladePosNed(durTimeIdx,2), -estimatedBladePosNed(durTimeIdx,3), 'b.', 'LineWidth', 1.0); grid on; hold on;
    plot3(posCmdNed(durTimeIdx,1), posCmdNed(durTimeIdx,2), -posCmdNed(durTimeIdx,3), '-r', 'LineWidth', 1.0); grid on; hold on;
    plot3(posNed(durTimeIdx,1), posNed(durTimeIdx,2), -posNed(durTimeIdx,3), '-k', 'LineWidth', 1.5); grid on; hold on;

    xlabel("Position(N)");
    ylabel("Position(E)");
    zlabel("Position(Alt)");
    axis ij;
    axis equal;
    legend('posCmd','vehicle position');
    
    % trajectory following performance analysis
    posErrNed = posCmdNed(durTimeIdx,:) - posNed(durTimeIdx,:);
    posDotT = dot(posErrNed, trajUnitVectorT((durTimeIdx),:),2)';
    posDotN = dot(posErrNed, trajUnitVectorN((durTimeIdx),:),2)';
    posDotB = dot(posErrNed, trajUnitVectorB((durTimeIdx),:),2)';
    posDotTNB = [posDotT; posDotN; posDotB];
    posDotTNBMean = mean(posDotTNB,2)
    posDotTNBStdDev = std(posDotTNB,0,2)
    posDotTNBMax = max(abs(posDotTNB(:,:)),[],2)
    
    figure(25);
    plot(time, trajTimeCur); hold on; 
    plot(time, trajTimeMax); grid on;
    legend('trajTimeCur','trajTimeMax')
end


%%
% Situation analysis in 3D
if situationAnalysis3d_show
    
    % find map for flight time index to payload image index
    plImageList = dir(gdLogImageFolder + "/*.jpg");
    plJsonList = dir(gdLogImageFolder + "/*.json");
    hourTimeHHmmssSS = hourTime;
    hourTimeHHmmssSS.Format = 'HHmmssSS';
    timeIdx2PlImageIdx = zeros(size(hourTimeHHmmssSS,1), 1);
    hourTimeFeature = str2double(string(hourTimeHHmmssSS));
    
    for i=1:size(plImageList,1)
        imageListNameStruct = split(plImageList(i).name, '.');
        imageListName(i) = string(imageListNameStruct{1});
        imageListNameStruct = split(imageListName(i), '_');
        imageListFeature(i,1) = str2double(imageListNameStruct(2)+imageListNameStruct{3}(1:2));
    end
    
    for i=1:size(timeIdx2PlImageIdx,1)
        [res, idx] = min(abs(hourTimeFeature(i) - imageListFeature));
        timeIdx2PlImageIdx(i) = idx;
    end
    
    % draw base figures
    sa3d_fig = figure(15);
    set(sa3d_fig, 'Position', [100, 50, 1200, 700]);
    axis3dplot = axes('Parent', sa3d_fig, 'position', [0.05, 0.2, 0.5, 0.75]);
    hold(axis3dplot, 'on');
    axisImage = axes('Parent', sa3d_fig, 'position', [0.6, 0.55, 0.35, 0.4]);

    % allocate GUIs
    sa3d_slider = uicontrol('Parent', sa3d_fig, 'Style', 'slider', 'Position', [80,50,1000,20], ...
                  'value', 1, 'min', 1, 'max', size(timeSec, 1), ...
                  'SliderStep', [1/size(timeSec,1), 10/size(timeSec,1)] );
    sa3d_bgcolor = sa3d_fig.Color;
    sa3d_text1 = uicontrol('Parent', sa3d_fig, 'Style', 'text', 'Position', [45,43,30,25], ...
                    'String', string(min(timeSec)),'BackgroundColor', sa3d_bgcolor);
    sa3d_text2 = uicontrol('Parent', sa3d_fig, 'Style', 'text', 'Position', [1080,43,70,25], ...
                    'String', string(max(timeSec)),'BackgroundColor', sa3d_bgcolor);
    sa3d_text3 = uicontrol('Parent', sa3d_fig, 'Style', 'text', 'Position', [550,25,100,25], ...
                    'String', '0.000', 'BackgroundColor', sa3d_bgcolor);

    % prepare data
    sa3dAxisHandle.plot3d = axis3dplot;
    sa3dAxisHandle.image = axisImage;
    sa3dGuiObject.textTime = sa3d_text3;
    sa3dData.folderName = gdLogImageFolder;
    sa3dData.time = time;
    sa3dData.fTimeLocal_Sec = fTimeLocal_Sec;
    sa3dData.vPosNed = posNed;
    sa3dData.vRpy = rpy;
    sa3dData.bdLidarDist = bdLidarDist;
    sa3dData.bdLidarAngle = bdLidarAngle;
    sa3dData.gimbalRpy = gimbalRPY;
    sa3dData.posCmdNed = posCmdNed;
    sa3dData.timeIdx2PlImageIdx = timeIdx2PlImageIdx;
    sa3dData.plImageList = plImageList;
    sa3dData.trajUvT = trajUnitVectorT;
    sa3dData.trajUvN = trajUnitVectorN;
    sa3dData.trajUvB = trajUnitVectorB;
    sa3dData.trajCmdTnb = trajCmdTnb;
    sa3dData.timeFormat = timeFormat;
    sa3dData.plJsonList = plJsonList;
    
    sa3d_slider.Callback = {@situationAnalysis3d, sa3dAxisHandle, sa3dGuiObject, sa3dData};
end


%%
% Remote Controller Input
if remoteController_show
    figure(16); 
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("Remote controller input");
    
    subFigRc(1) = subplot(5,1,1);
    plot(time, jobSeq, '-k', 'LineWidth', 1.5); hold on; 
    plot(time, fcMcMode, '--b', 'LineWidth', 1.5); hold on;
    plot(time, ctrlSetpointType, '-.r', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    ylim([0 22]);
    legend("jobSeq","fcMcMode","ctrlSetpointType");
    ylabel("Index");
    grid on;

    subFigRc(2) = subplot(5,1,2);
    plot(time, rcPitch, '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    legend("rcPitch");
    ylabel("rcPitch (%)");
    grid on;

    subFigRc(3) = subplot(5,1,3);
    plot(time, rcRoll, '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    legend("rcRoll");
    ylabel("rcRoll (%)");
    grid on;

    subFigRc(4) = subplot(5,1,4);
    plot(time, rcYaw, '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    legend("rcYaw");
    ylabel("rcYaw (%)");
    grid on;
    
    subFigRc(5) = subplot(5,1,5);
    plot(time, rcThrottle, '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    legend("rcThrottle");
    ylabel("rcThrottle (%)");
    grid on;

    linkaxes([subFigRc(1), subFigRc(2), subFigRc(3), subFigRc(4), subFigRc(5)], 'x');
    xlim([sTime eTime]);
end

%%
% Gimbal Angle
if gimbalAngle_show
    figure(17); 
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("Gimbal angle response");
    
    subFigGimbal(1) = subplot(4,1,1);
    plot(time, jobSeq, '-k', 'LineWidth', 1.5); hold on; 
    plot(time, fcMcMode, '--b', 'LineWidth', 1.5); hold on;
    plot(time, ctrlSetpointType, '-.r', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
    ylim([0 22]);
    legend("jobSeq","fcMcMode","ctrlSetpointType");
    ylabel("Index");
    grid on;

    subFigGimbal(2) = subplot(4,1,2);
    plot(time, gimbalRPY(:,1), '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    legend("gimbal roll");
    ylabel("gimbal roll (deg)");
    grid on;

    subFigGimbal(3) = subplot(4,1,3);
    plot(time, gimbalRPY(:,2), '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    legend("gimbal pitch");
    ylabel("gimbal pitch (deg)");
    grid on;

    subFigGimbal(4) = subplot(4,1,4);
    plot(time, gimbalRPY(:,3), '-k', 'LineWidth', 1.2);
    xtickformat(strTimeXtickFormat);
    legend("gimbal yaw");
    ylabel("gimbal yaw (deg)");
    grid on;

    linkaxes([subFigGimbal(1), subFigGimbal(2), subFigGimbal(3), subFigGimbal(4)], 'x');
    xlim([sTime eTime]);
end

%%
% bladeTravelDistance
if bladeTravelDistance_show
    figure(18); 
    set(gcf, 'Position', [100, -20, 800, 800]);
    suptitle("Blade Travel Distance");
    
    subFigGimbal(1) = subplot(1,1,1);
    plot(time, bladeTravelDistance, '-k', 'LineWidth', 1.5); hold on; 
%     plot(time, fcMcMode, '--b', 'LineWidth', 1.5); hold on;
%     plot(time, ctrlSetpointType, '-.r', 'LineWidth', 1.5);
    xtickformat(strTimeXtickFormat);
%     ylim([0 22]);
    legend("bladeTravelDistance");
    ylabel("Distance");
    grid on;
    xlim([sTime eTime]);
end

%% Mission Time Analysis
if missionTime_show
    
    Idiff = @(x) [1, find(diff(x)-1)+1, length(x)];
    
%     JobSeqArr = [1, 4, 6, 9, 10]; % Blade inspection job for 5-1 side mission
%     JobSeqArr = [1, 4, 6, 10, 11]; % Blade inspection job for 5-2 side mission

%     JobSeqArr = [1, 4, 6, 9, 10]; % Blade inspection job for 4-1 side mission
%     JobSeqArr = [1, 4, 6, 9, 10]; % Blade inspection job for 4-2 side mission

%     TotalJob = (1:16); % Total job for 6_1 pitch mission
%     TotalJob = (1:17); % Total job for 6_2 pitch, 6_1 no pitch mission
    TotalJob = (1:20); % Total job for 6_2 no pitch mission

%   6_1 no pitch
%     JobSeqArr = [1, 4, 6, 9, 14, 17];  % Blade inspection job for 6_1 no pitch mission
%     CamTranSeqArr = [5, 10, 13]; % Cam Transition job for 6_1 no pitch mission

%   6_2 no pitch
    JobSeqArr = [1, 5, 10, 14, 16 ,19];  % Blade inspection job for 6_2 no pitch mission
    CamTranSeqArr = [6, 9, 15]; % Cam Transition job for 6_2 no pitch mission

%   6_1 / 6_2 pitch
%     JobSeqArr = [1, 4, 6, 10, 12, 15];  % Blade inspection job for 6_1 / 6_2 pitch mission
%     CamTranSeqArr = [5, 11]; % Cam Transition job for 6_1 / 6_2 pitch mission

    WpJobSeqArr = setdiff(setdiff(TotalJob, JobSeqArr), CamTranSeqArr);      
    WpIdx = Idiff(WpJobSeqArr); 
    WpLdx = diff(WpIdx);
    WpLdx(end) = WpLdx(end)+1;
    WpIdx = WpIdx(1:end-1);
    WpJobTime = zeros(1, length(WpJobSeqArr));
    
    Job.stdJobTime = zeros(1, length(JobSeqArr));
    Job.camTranJobTime = zeros(1, length(CamTranSeqArr));
    Job.wpJobTime = zeros(1, length(WpIdx));
    Job.wpEachJobTime = zeros(1, length(WpJobSeqArr));
    Job.stdJobNumber = length(JobSeqArr);
    Job.camTranJobNumber = length(CamTranSeqArr);
    Job.wpJobNumber = length(WpJobSeqArr);
    
    
    for i = 1:size(JobSeqArr,2)
        JobNumb = JobSeqArr(i);

        durJob = find(jobSeq == JobNumb);    
        durJob = find(fcMcMode(min(durJob):max(durJob)) == 2);
        Job.stdJobTime(i) = length(durJob)*0.02;
    end
    
    for i = 1:size(CamTranSeqArr,2)
        JobNumb = CamTranSeqArr(i);

        durJob = find(jobSeq == JobNumb);    
        Job.camTranJobTime(i) = length(durJob)*0.02;
    end
    
    for i = 1:size(WpJobSeqArr,2)
        JobNumb = WpJobSeqArr(i);
        
        durJob = find(jobSeq == JobNumb);
        durJob = find(fcMcMode(min(durJob):max(durJob)) == 2);
        Job.wpEachJobTime(i) = length(durJob)*0.02;
    end
    
    for i = 1: length(WpIdx)
        Job.wpJobTime(i) = sum(Job.wpEachJobTime(WpIdx(i):WpIdx(i)+WpLdx(i)-1));
    end
    Job.totalStdJobTime = sum(Job.stdJobTime);
    Job.totalCamTranJobTime = sum(Job.camTranJobTime);
    Job.totalWpJobTime = sum(Job.wpJobTime);
    Job.totalTime = sum([Job.totalStdJobTime, Job.totalCamTranJobTime, Job.totalWpJobTime]);
    
end

%% 
% Analysis for each blade inspection job for 5-1, 5-2 missions
if analEachJob_show
    
%     JobSeqArr = [1, 4, 6, 9, 10]; % Blade inspection job for 5-1 side mission
    % JobSeqArr = [1, 4, 6, 10, 11]; % Blade inspection job for 5-2 side mission
    JobSeqArr = [1, 4, 6, 9, 10]; % Blade inspection job for 4-1 side mission
%     JobSeqArr = [1, 4, 6, 9, 10]; % Blade inspection job for 4-2 side mission
    % JobSeqArr = [1]; % test

    ctrlSpUvw_t = ctrlSp;
    ctrlSpType_t = ctrlSetpointType;
    for i = 1:size(JobSeqArr,2)
        JobNumb = JobSeqArr(i);

        durJob = find(jobSeq == JobNumb);

        ctrlSpType_t([min(durJob):max(durJob)], 1);

        durMinIdx = (min(durJob) + 10);
    %     durMax = min(durJob) + min(find(ctrlSpType_t([min(durJob):max(durJob)]) ~= 8)) - 10;
        durMaxIdx = (min(durJob) + min(find(fcMcMode([(min(durJob)+10):max(durJob)]) == 1)) - 10);

        durIdx = durMinIdx:durMaxIdx;

        errVelNav = velNav(durIdx,:) - velCmdNav(durIdx,:);
        rmse_errVelNav(:,:,i) = sqrt(mean(errVelNav.^2));

        figure; set(gcf, 'Position', [100, -20, 800, 800]);
        supTitle_array = "Velocity Responses (UVW) for Job # " + JobSeqArr(i);
        suptitle(supTitle_array);

        subFigVel(1) = subplot(3,1,1);
        plot(time(durIdx), velCmdNav(durIdx,1), '--r', 'LineWidth', 1.5); hold on;
        plot(time(durIdx), velNav(durIdx,1), '-k', 'LineWidth', 1.2); 
        xtickformat(strTimeXtickFormat);
        xlim([time(durMinIdx) time(durMaxIdx)]);
        ylim([-1.5 1.5]);
        legend("Velocity Sp (NavX)", "Velocity (NavX)");
        ylabel("Velocity (Nav X) (m/s)");
        grid on;

        subFigVel(2) = subplot(3,1,2);
        plot(time(durIdx), velCmdNav(durIdx,2), '--r', 'LineWidth', 1.5); hold on;
        plot(time(durIdx), velNav(durIdx,2), '-k', 'LineWidth', 1.2);
        xtickformat(strTimeXtickFormat);
        xlim([time(durMinIdx) time(durMaxIdx)]);
        ylim([-1.5 1.5]);
        legend("Velocity Sp (NavY)", "Velocity (NavY)");
        ylabel("Velocity (Nav Y) (m/s)");
        grid on;

        subFigVel(3) = subplot(3,1,3);

        plot(time(durIdx), velCmdNav(durIdx,3), '--r', 'LineWidth', 1.5); hold on;
        plot(time(durIdx), velNav(durIdx,3), '-k', 'LineWidth', 1.2);
        xtickformat(strTimeXtickFormat);
        xlim([time(durMinIdx) time(durMaxIdx)]);
        ylim([-1.5 1.5]);
        legend("Velocity Sp (NavZ)", "Velocity (NavZ)");
        ylabel("Velocity (Nav Z) (m/s)");
        xlabel("time (sec)");
        grid on;

        linkaxes([subFigVel(1), subFigVel(2), subFigVel(3)], 'x');
    end
    rmse_errVelNav
    
end

