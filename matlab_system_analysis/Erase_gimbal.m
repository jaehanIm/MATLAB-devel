% close all

% data = readtable('/home/jaehan/log/210723/gdLog_210723_114738.csv');
data = readtable('/home/jaehan/log/210802/gdLog_210805_102755.csv');

time = data.rosTime;
time = seconds(time);
time = time - time(1);
fcMode = data.fcMcMode;

ns = 1; nf = 1; nsm = 1; nfm = 1;
for i = 2:length(fcMode)
    if fcMode(i) ~= fcMode(i-1) && fcMode(i) == 2
        testStartFlag(ns) = i; ns = ns+1;
    elseif fcMode(i) ~= fcMode(i-1) && fcMode(i-1) == 2
        testFinishFlag(nf) = i; nf = nf+1;
    end
    if fcMode(i-1) == 0 && fcMode(i) == 1
        missionStartFlag(nsm) = i; nsm = nsm + 1;
    elseif fcMode(i-1) == 1 && fcMode(i) == 0
        missionFinishFlag(nfm) = i; nfm = nfm + 1;
    end
end

data(testStartFlag(8):testFinishFlag(8)-1,:) = [];

time = data.rosTime;
time = seconds(time);
time = time - time(1);
fcMode = data.fcMcMode;

ns = 1; nf = 1; nsm = 1; nfm = 1;
testStartFlag = [];
testFinishFlag = [];
missionStartFlag = [];
missionFinishFlag = [];
for i = 2:length(fcMode)
    if fcMode(i) ~= fcMode(i-1) && fcMode(i) == 2
        testStartFlag(ns) = i; ns = ns+1;
    elseif fcMode(i) ~= fcMode(i-1) && fcMode(i-1) == 2
        testFinishFlag(nf) = i; nf = nf+1;
    end
    if fcMode(i-1) == 0 && fcMode(i) == 1
        missionStartFlag(nsm) = i; nsm = nsm + 1;
    elseif fcMode(i-1) == 1 && fcMode(i) == 0
        missionFinishFlag(nfm) = i; nfm = nfm + 1;
    end
end





%%

figure(1)
clf
subplot(2,1,1)
title('Gimbal Roll Response')
hold on
plot(time,data.gimbalRpy_deg_0)
plot(time,data.gimbalRpyCmd_deg_0)
% plot(time(2:end),diff(data.rosTime)*100)
plot(time,data.fcMcMode)
grid on
subplot(2,1,2)
title('Gimbal Pitch Response')
hold on
plot(time,data.gimbalRpy_deg_1)
plot(time,data.gimbalRpyCmd_deg_1)
% plot(time(2:end),diff(data.rosTime))
plot(time,data.fcMcMode)
grid on


%% TF estimation
numStepRoll = 9;
numStepPitch = 9;
tfSetRoll = {};
tfSetPitch = {};
tfSetSweepRoll = {};
tfSetSweepPitch = {};

for i = 2:numStepRoll+1
    cmd_ = data.gimbalRpyCmd_deg_0(testStartFlag(i):testFinishFlag(i)-1);
    resp_ = data.gimbalRpy_deg_0(testStartFlag(i):testFinishFlag(i)-1);
    tfSetRoll{i-1} = tfestimator(cmd_, resp_);
end

for i = 1:numStepPitch
    cmd_ = data.gimbalRpyCmd_deg_1(testStartFlag(i+10):testFinishFlag(i+10)-1);
    resp_ = data.gimbalRpy_deg_1(testStartFlag(i+10):testFinishFlag(i+10)-1);
    tfSetPitch{i} = tfestimator(cmd_, resp_);
end

cmd_ = data.gimbalRpyCmd_deg_0(testStartFlag(20):testFinishFlag(20)-1);
resp_ = data.gimbalRpy_deg_0(testStartFlag(20):testFinishFlag(20)-1);
tfSetSweepRoll = tfestimator(cmd_, resp_);

cmd_ = data.gimbalRpyCmd_deg_1(testStartFlag(21):testFinishFlag(21)-1);
resp_ = data.gimbalRpy_deg_1(testStartFlag(21):testFinishFlag(21)-1);
tfSetSweepPitch = tfestimator(cmd_, resp_);

%% RiseTime plot
rollRT = [];
pitchRT = [];
for i = 1:numStepRoll
    rollRT(i) = stepinfo(tf(tfSetRoll{i})).RiseTime;
    pitchRT(i) = stepinfo(tf(tfSetPitch{i})).RiseTime;
end

figure(2)
clf
hold on
grid on
plot(rollRT(3:end))
plot(pitchRT(3:end))

%% Step plot
figure(3)
clf
for i = 1:9
    subplot(3,3,i)
    step(tfSetRoll{i})
    grid on
    titleText = ['testNum : ',num2str(i)];
    title(titleText)
end

figure(4)
clf
for i = 1:9
    subplot(3,3,i)
    step(tfSetPitch{i})
    grid on
    titleText = ['testNum : ',num2str(i)];
    title(titleText)
end

figure(5)
clf
for i = 1:9
    subplot(3,3,i)
    margin(tfSetRoll{i})
    grid on
    titleText = ['testNum : ',num2str(i)];
    title(titleText)
end

figure(6)
clf
for i = 1:9
    subplot(3,3,i)
    margin(tfSetPitch{i})
    grid on
    titleText = ['testNum : ',num2str(i)];
    title(titleText)
end

figure(7)
clf
subplot(2,1,1)
step(tfSetSweepRoll)
subplot(2,1,2)
margin(tfSetSweepRoll)

figure(8)
clf
subplot(2,1,1)
step(tfSetSweepPitch)
subplot(2,1,2)
margin(tfSetSweepPitch)

%%
function tf = tfestimator(input, output)
timeSeries = iddata(output,input,0.02);
sysR = tfest(timeSeries, 3, 2, nan);
tf = sysR;
end