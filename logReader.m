function [time,out,acc,gimbal] = logReader(gdLogFileName)
d2r = pi/180;

useMinuteFormat = true;
useHourFormat = false;
gdLog = importfile(gdLogFileName);

%% Loading data
gdLog = importfile(gdLogFileName);
data = gdLog;
timeValidIdx = max(find(gdLog.rosTime < 1514764800)); % 2018/01/01 == 1514764800
data((1:timeValidIdx),:) = [];
dataSize = size(data.rosTime);
time = data.rosTime;
time = time(:,1) - time(1,1);
time = seconds(time);
dateStart = datetime(data.rosTime(2),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
if(useMinuteFormat)
    strTimeXtickFormat = 'mm:ss.SSS';
else
    strTimeXtickFormat = 'ss.SSS';
end

if(useHourFormat)
    time = time + dateStart;
    strTimeXtickFormat = 'HH:mm:ss';
end

jobSeq = data.jobSeq;
rpy(:,1) = data.rpy_0;
rpy(:,2) = data.rpy_1;
rpy(:,3) = data.rpy_2;
pos(:,1) = data.posNed_0;
pos(:,2) = data.posNed_1;
pos(:,3) = data.posNed_2;
rpyRad = rpy*d2r;
fcMcMode = data.fcMcMode;
acc(:,1) = data.accBody_0;
acc(:,2) = data.accBody_1;
acc(:,3) = data.accBody_2;
vel(:,1) = data.velNed_0;
vel(:,2) = data.velNed_1;
vel(:,3) = data.velNed_2;
gimbal(:,1) = data.gimbalRPY_0;
gimbal(:,2) = data.gimbalRPY_1;
gimbal(:,3) = data.gimbalRPY_2;

out = [fcMcMode,jobSeq,rpy,pos,vel];

sTime = min(time);
eTime = max(time);
timeDurTotal = [sTime eTime];

%% Main.m
end
