addpath('./..')

%% Load data
% 221110_1 - C #1
% gdLog = readtable('/home/jaehan/log/221110_113922/gdLog_221110_113922.csv');
% imuData = readtable('/home/jaehan/log/221110_113922/aSensorImu_221110_113922.csv');

% 221110_2 - C #2
% gdLog = readtable('/home/jaehan/log/221110_121923/gdLog_221110_121923.csv');
% imuData = readtable('/home/jaehan/log/221110_121923/aSensorImu_221110_121923.csv');

% 221110_3 - body 1_!
% gdLog = readtable('/home/jaehan/log/221110_123033/gdLog_221110_123033.csv');
% imuData = readtable('/home/jaehan/log/221110_123033/aSensorImu_221110_123033.csv');

% 221115_4 - D
% gdLog = readtable('/home/jaehan/log/221115_111538/gdLog_221115_111538.csv');
% imuData = readtable('/home/jaehan/log/221115_111538/aSensorImu_221115_111538.csv');

% 221115_5 - body 1_2
% gdLog = readtable('/home/jaehan/log/221115_112930/gdLog_221115_112930.csv');
% imuData = readtable('/home/jaehan/log/221115_112930/aSensorImu_221115_112930.csv');

% 221201_1 - nowind
% gdLog = readtable('/home/jaehan/log/221201_111806/gdLog_221201_111806.csv');
% imuData = readtable('/home/jaehan/log/221201_111806/aSensorImuLog_221201_111806.csv');

% 221201_1 - midwind
% gdLog = readtable('/home/jaehan/log/221201_115426/gdLog_221201_115426.csv');
% gdLog(1:13516,:) = [];
% imuData = readtable('/home/jaehan/log/221201_115426/aSensorImuLog_221201_115426.csv');

% 221201_1 - highwind
% gdLog = readtable('/home/jaehan/log/221201_120533/gdLog_221201_120533.csv');
% imuData = readtable('/home/jaehan/log/221201_120533/aSensorImuLog_221201_120533.csv');

% 230212 - wind
gdLog = readtable('/data/230228_121753/gdLog_230228_121753.csv');
imuData = readtable('/data/230228_121753/aSensorImuLog_230228_121753.csv');

jobInterest = 5;
freqL = [0 5  30  75 150 800 1200];
freqH = [5 30 75 150 800 1200 2000];

for i = 1:size(imuData,1)
    if imuData.rosTime(i) > 100000
        thres = i;
        break;
    end
end
imuData(1:thres-1,:) = [];

% Job index parser
jobStartIdx = [];
jobEndIdx = [];
landIdx = [];
for i = 2:size(gdLog,1)
    if gdLog.fcMcMode(i) == 2 && gdLog.fcMcMode(i-1) == 1
        jobStartIdx = horzcat(jobStartIdx,i);
    elseif (gdLog.fcMcMode(i) == 1 && gdLog.fcMcMode(i-1) == 2) ||(gdLog.fcMcMode(i) == 0 && gdLog.fcMcMode(i-1) == 2)
        jobEndIdx = horzcat(jobEndIdx,i);
    elseif gdLog.fcMcMode(i) == 255 && gdLog.fcMcMode(i-1) == 0
        landIdx = i;
    end
end

global start finish missStartTime missEndTime imuTimeS

% time synchronization
imuTime = datetime(imuData.rosTime,'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
gdTime = datetime(gdLog.rosTime(1:end),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
imuTimeDelay = seconds(imuTime(1) - gdTime(1));
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);
imuTimeS = seconds(imuTime - imuTime(1)) + imuTimeDelay;
% imuTimeS = seconds(imuTime - imuTime(1)) + 51.9 + 1.67;

% Assign data - !!! sensor attitude aware !!!
% acc_0 = imuData.acc_mpss_0;
% acc_1 = imuData.acc_mpss_1;
% acc_2 = imuData.acc_mpss_2;
% gyro_0 = imuData.gyro_dps_0;
% gyro_1 = imuData.gyro_dps_1;
% gyro_2 = imuData.gyro_dps_2;

acc_0 = imuData.acc_mpss_2;
acc_1 = imuData.acc_mpss_1;
acc_2 = imuData.acc_mpss_0;
gyro_0 = imuData.gyro_dps_2;
gyro_1 = imuData.gyro_dps_1;
gyro_2 = imuData.gyro_dps_0;
% 
% acc_0 = -imuData.acc_mpss_2;
% acc_1 = -imuData.acc_mpss_1;
% acc_2 = imuData.acc_mpss_0;
% gyro_0 = -imuData.gyro_dps_2;
% gyro_1 = -imuData.gyro_dps_1;
% gyro_2 = imuData.gyro_dps_0;

Fs = round(1/mean(diff(imuTimeS)));
L = size(imuData,1);

missIdxTemp = find(diff(gdLog.fcMcMode)~=0);
missStartTime = gdTimeS([missIdxTemp(2), missIdxTemp(4), missIdxTemp(6), missIdxTemp(8), missIdxTemp(10)]);
missEndTime = gdTimeS([missIdxTemp(3), missIdxTemp(5), missIdxTemp(7), missIdxTemp(9), missIdxTemp(11)]);

setImuIdx(jobInterest)

disp("Data Loading Complete!")

%% pre processing
acc_0 = detrend(acc_0,1);
acc_1 = detrend(acc_1,1);
acc_2 = detrend(acc_2,2);
gyro_0 = detrend(gyro_0,1);
gyro_1 = detrend(gyro_1,1);
gyro_2 = detrend(gyro_2,1);

%% Raw plot
figure(2)
clf

subplot(6,1,1)
plot(imuTimeS(start:finish),gyro_0(start:finish))
title('gyro R')
xlabel('time[s]');
ylabel('[rad/s]');
hold on
ylim([-10*std(gyro_0) 10*std(gyro_0)])
xlim([imuTimeS(start) imuTimeS(finish)])

subplot(6,1,2)
plot(imuTimeS(start:finish),gyro_1(start:finish))
title('gyro P')
xlabel('time[s]')
ylabel('[rad/s]');
hold on
ylim([-10*std(gyro_1) 10*std(gyro_1)])
xlim([imuTimeS(start) imuTimeS(finish)])

subplot(6,1,3)
plot(imuTimeS(start:finish),gyro_2(start:finish))
title('gyro Y')
xlabel('time[s]')
ylabel('[rad/s]');
hold on
ylim([-10*std(gyro_2) 10*std(gyro_2)])
xlim([imuTimeS(start) imuTimeS(finish)])

subplot(6,1,4)
plot(imuTimeS(start:finish),acc_0(start:finish))
title('acc X')
xlabel('time[s]');
ylabel('[m/s^2]');
hold on
ylim([-10*std(acc_0) 10*std(acc_0)])
xlim([imuTimeS(start) imuTimeS(finish)])

subplot(6,1,5)
plot(imuTimeS(start:finish),acc_1(start:finish))
title('acc Y')
xlabel('time[s]')
ylabel('[m/s^2]');
hold on
ylim([-10*std(acc_1) 10*std(acc_1)])
xlim([imuTimeS(start) imuTimeS(finish)])

subplot(6,1,6)
plot(imuTimeS(start:finish),acc_2(start:finish))
title('acc Z')
xlabel('time[s]')
ylabel('[m/s^2]');
hold on
ylim([-10*std(acc_2) 10*std(acc_2)])
xlim([imuTimeS(start) imuTimeS(finish)])

%% FDI

t = imuTime;
ax = acc_0;
ay = acc_1;
az = acc_2;
L = length(ax);

% regen engine
u = FDI(ax,Fs);
v = FDI(ay,Fs);
w = FDI(az,Fs);
    
% HPF
filterFreqUpper = 3; %hz
U = fft(u);
V = fft(v);
W = fft(w);

df = L/Fs;
kernel = ones(L,1);
kernel(1:ceil(filterFreqUpper*df)) = 0;
kernel(L-ceil(filterFreqUpper*df)+1:end) = 0;

U_filter = U .* kernel;
V_filter = V .* kernel;
W_filter = W .* kernel;

u = real(ifft(U_filter));
v = real(ifft(V_filter));
w = real(ifft(W_filter));

figure(4)
clf
subplot(2,1,1)
hold on
plot(imuTimeS,v)
grid on
xlabel('time[s]')
ylabel('v [m/s]');
plot(gdTimeS,gdLog.velUvw_mps_1)
title('velocity trend comparison')
legend('integrated v (HP filtered)','gdLog v')


%% Partial FDI
for i = 1:length(gdTimeS)-1
    interest = find(imuTimeS>=gdTimeS(i) & imuTimeS<=gdTimeS(i+1));
    
    gimbalR = gdLog.gimbalRpy_deg_0(i);
    gimbalP = gdLog.gimbalRpy_deg_1(i);
    gimbalY = 0;
    dcmI2G = angle2dcm(gimbalY, gimbalR, gimbalP,'ZXY');
    
    timeF = gdTimeS(i);
    timeB = gdTimeS(i+1);
    gduF = gdLog.velUvw_mps_0(i);
    gduB = gdLog.velUvw_mps_0(i+1);
    gdvF = gdLog.velUvw_mps_1(i);
    gdvB = gdLog.velUvw_mps_1(i+1);
    gdwF = gdLog.velUvw_mps_2(i);
    gdwB = gdLog.velUvw_mps_2(i+1);
    
    gdVelF = dcmI2G * [gduF;gdvF;gdwF];
    gdVelB = dcmI2G * [gduB;gdvB;gdwB];
    
    gduF = gdVelF(1); gdvF = gdVelF(2); gdwF = gdVelF(3);
    gduB = gdVelB(1); gdvB = gdVelB(2); gdwB = gdVelB(3);
    
    interestRatio = (imuTimeS(interest)-timeF)/(timeB-timeF);
%     max(interestRatio)
%     min(interestRatio)
%     disp(cut' + num2str(i))
    
    u(interest) = u(interest) + (1-interestRatio) * gduF + interestRatio * gduB;
    v(interest) = v(interest) + (1-interestRatio) * gdvF + interestRatio * gdvB;
    w(interest) = w(interest) + (1-interestRatio) * gdwF + interestRatio * gdwB;
end

figure(4)
subplot(2,1,2)
hold on
plot(imuTimeS,v)
grid on
xlabel('time[s]')
ylabel('v [m/s]');
plot(gdTimeS,gdLog.velUvw_mps_1)
title('velocity trend comparison')
legend('integrated v (HP filtered)','gdLog v')

disp('Partial FDI complete!')

figure(55)
subplot(3,1,1)
hold on
plot(imuTimeS,u)
plot(gdTimeS,gdLog.velUvw_mps_0)
subplot(3,1,2)
hold on
plot(imuTimeS,v)
plot(gdTimeS,gdLog.velUvw_mps_1)
subplot(3,1,3)
hold on
plot(imuTimeS,w)
plot(gdTimeS,gdLog.velUvw_mps_2)
%% data selection
% y = gyro_1;

%% FFT plot
% figure(2)
% title('FFT')
% clf
% [freq, fftResult, psdResult] = data2fftpsd(u,Fs);
% % plot(freq,fftResult);
% plot(freq,fftResult);
% hold on
% [freq, fftResult, psdResult] = data2fftpsd(v,Fs);
% % plot(freq,fftResult);
% plot(freq,fftResult);
% grid on
% xlabel('hz');
% ylabel('magnitude')
% legend('u','v')
% xlim([0.1 2000])
% 
% figure(3)
% title('PSD')
% clf
% [freq, fftResult, psdResult] = data2fftpsd(gyro_0,Fs);
% plot(freq,fftResult);
% hold on
% [freq, fftResult, psdResult] = data2fftpsd(gyro_1,Fs);
% plot(freq,fftResult);
% grid on
% xlabel('hz');
% ylabel('psd')
% legend('gyro_R','gyro_P')
% xlim([0.1 2000])

%% Spectrogram
% y = acc_1;
y = acc_1;
timeStep = 1; %[s]

dspec = timeStep * Fs;
specNum = floor(L/dspec);
% spectrogram = zeros(specNum,L);
spectrogram = [];

for i = 1:specNum
    y_spec = y(dspec*(i-1)+1:dspec*i);
    [freq,fftResult,~] = data2fftpsd(y_spec,Fs);
    spectrogram(i,1:length(freq)) = fftResult';
end

spectrogram(:,length(freq)+1:end) = [];

spectrogram(1:32,:) = [];
% spectrogram = log(spectrogram);

figure(8)
clf
imagesc([freq(1) freq(end)],[imuTimeS(1) imuTimeS(end)],spectrogram(1:end,1:end));
% imagesc([freq(1) freq(end)],[missStartTime(1) missEndTime(end)],spectrogram(1:end,1:end));
xlabel('freq [hz]')
ylabel('time [s]')
xlim([0 400])
title('Spectrogram Analysis [dt = 1.0s]')
colorbar
colormap("turbo")

figure(88)
clf
sgtitle('spectrogram')
colormap("turbo")
subplot(1,4,1)
imagesc([0 30],[imuTimeS(1) imuTimeS(end)],spectrogram(:,0*timeStep+1:30*timeStep));
colorbar
xlabel('hz')
ylabel('time [s]')
title('0~30hz')
hold on


subplot(1,4,2)
imagesc([30 75],[imuTimeS(1) imuTimeS(end)],spectrogram(:,30*timeStep:75*timeStep));
hold on
xlabel('hz')
ylabel('time [s]')
title('30~75hz')
colorbar

subplot(1,4,3)
imagesc([75 150],[imuTimeS(1) imuTimeS(end)],spectrogram(:,75*timeStep:150*timeStep));
hold on
xlabel('hz')
ylabel('time [s]')
title('75~150')
colorbar

subplot(1,4,4)
imagesc([150 400],[imuTimeS(1) imuTimeS(end)],spectrogram(:,150*timeStep:400*timeStep));
hold on
xlabel('hz')
ylabel('time [s]')
title('150~400hz')
colorbar

figure(888)
clf
subplot(3,1,1)
hold on
grid on
imagesc([imuTimeS(1) imuTimeS(end)],[0 200],spectrogram(:,0*timeStep+1:200*timeStep)');
colorbar
xlabel('time [s]')
ylabel('hz')
title('0~200hz')
hold on
xlim([imuTimeS(1) imuTimeS(end)])
ylim([0 200])

subplot(3,1,2)
hold on
grid on
% plot(gdTimeS,gdLog.pqr_dps_0)
% plot(gdTimeS,gdLog.pqr_dps_1)
% plot(gdTimeS,gdLog.pqr_dps_2)
% plot(gdTimeS,gdLog.pqrBody_dps_0)
% plot(gdTimeS,gdLog.pqrBody_dps_1)
% plot(gdTimeS,gdLog.pqrBody_dps_2)
xlim([imuTimeS(1) imuTimeS(end)])
colorbar


disp("Spectrogram plotting complete!")


%% Vibration component analysis - pqr
p = gyro_0;
q = gyro_1;
r = gyro_2;

P = fft(p);
Q = fft(q);
R = fft(r);

df = L/Fs;

filteredFreq_p = [];
filteredFreq_q = [];
filteredFreq_r = [];

for i = 1:length(freqL)
    P_filter = P;
    Q_filter = Q;
    R_filter = R;
    
    kernel = zeros(L,1);
    kernel(max(1,freqL(i)*df):freqH(i)*df) = 1;
    kernel(L-freqH(i)*df:L-freqL(i)*df) = 1;
    
    P_filter = P_filter .* kernel;
    Q_filter = Q_filter .* kernel;
    R_filter = R_filter .* kernel;
    
    fftRegen_p = ifft(P_filter);
    fftRegen_q = ifft(Q_filter);
    fftRegen_r = ifft(R_filter);
    
    filteredFreq_p(i,:) = real(fftRegen_p);
    filteredFreq_q(i,:) = real(fftRegen_q);
    filteredFreq_r(i,:) = real(fftRegen_r);
end

vel = [];
x_comp = [];
y_comp = [];
for i = 1:length(freqL)
    x_comp(i,:) = filteredFreq_r(i,:) - 3 * filteredFreq_p(i,:) * cos(pi/4) / 9;
    y_comp(i,:) = filteredFreq_q(i,:) + 3 * filteredFreq_p(i,:) * sin(pi/4) / 9;
    vel(i,:) = sqrt(x_comp(i,:).^2 + y_comp(i,:).^2) * 180 / pi;
end

% totalVelocity = sum(vel,1);
% totalVelocity = sqrt(sum(filteredFreq_q).^2 + sum(filteredFreq_r).^2) * 180 / pi;
totalVelocity = sqrt(sum(x_comp).^2 + sum(y_comp).^2) * 180 / pi;
threshold = 24.8; % deg/s
residue = threshold - totalVelocity;
velBar = vertcat(vel,residue);

%% Vibration component analysis - uvw

U = fft(u);
V = fft(v);
W = fft(w);

df = L/Fs;

filteredFreq_u = [];
filteredFreq_v = [];
filteredFreq_w = [];

for i = 1:length(freqL)
    U_filter = U;
    V_filter = V;
    W_filter = W;
    
    kernel = zeros(L,1);
    kernel(max(1,freqL(i)*df):freqH(i)*df) = 1;
    kernel(L-freqH(i)*df:L-freqL(i)*df) = 1;
    
    U_filter = U_filter .* kernel;
    V_filter = V_filter .* kernel;
    W_filter = W_filter .* kernel;
    
    fftRegen_u = ifft(U_filter);
    fftRegen_v = ifft(V_filter);
    fftRegen_w = ifft(W_filter);
    
    filteredFreq_u(i,:) = real(fftRegen_u);
    filteredFreq_v(i,:) = real(fftRegen_v);
    filteredFreq_w(i,:) = real(fftRegen_w);
end

transvel = [];
for i = 1:length(freqL)
    transvel(i,:) = sqrt(filteredFreq_v(i,:).^2 + filteredFreq_w(i,:).^2);
end

% totalVelocity = sum(vel,1);
totalVelocityTrans = sqrt(sum(filteredFreq_v).^2 + sum(filteredFreq_w).^2);
thresholdTrans = 3.8; % m/s
residueTrans = thresholdTrans - totalVelocityTrans;
velBarTrans = vertcat(transvel,residueTrans);

%% vibration component plot angular

figure(9)
clf
subplot(2,1,1)
hold on
grid on
plot(imuTimeS(start:finish),totalVelocity(start:finish),'k','LineWidth',1)
plot([imuTimeS(start) imuTimeS(finish)],[24.8 24.8], 'r--','LineWidth',2)
xlim([imuTimeS(start) imuTimeS(finish)])
title('Pointing error velocity - angular component');
legend('pointing velocity','threshold');
xlabel('time[s]');
ylabel('deg/s');

subplot(2,1,2)
a = area(imuTimeS(start:finish),velBar(1:end-1,start:finish)');
xlim([imuTimeS(start) imuTimeS(finish)])
xlabel('time[s]')
ylabel('deg/s')
title('Contribution ratio area plot - angular component')
hold on
plot([imuTimeS(start) imuTimeS(finish)],[24.8 24.8], 'r:','LineWidth',2.4)
legend('DC-5hz','5-30hz','30-75hz','75-150hz','150-800hz','800-1200hz','1200> hz','threshold')
grid on

dpsViolationRate = size(find(totalVelocity >= threshold),2)/L * 100
dpsViolationRate = size(find(totalVelocity(start:finish) >= threshold),2)/(finish-start) * 100



count = 0;
for i= 1:L-20
if(~isempty(find(totalVelocity(i:i+20) > threshold)))
count = count + 1;
end
end
count/(L-20) * 100

%% vibration component plot translational

figure(10)
clf
subplot(2,1,1)
hold on
grid on
plot(imuTimeS(start:finish),totalVelocityTrans(start:finish),'k','LineWidth',1)
plot([imuTimeS(start) imuTimeS(finish)],[thresholdTrans thresholdTrans], 'r--','LineWidth',2)
xlim([imuTimeS(start) imuTimeS(finish)])
title('Pointing error velocity - translational component');
legend('pointing velocity','threshold');
xlabel('time[s]');
ylabel('m/s');
xlim([imuTimeS(start) imuTimeS(finish)])

subplot(2,1,2)

a = area(imuTimeS(start:finish),velBarTrans(1:end-1,start:finish)');
xlim([imuTimeS(start) imuTimeS(finish)])
xlabel('time[s]')
ylabel('m/s')
title('Contribution ratio area plot - translational component')
hold on
% a(8).FaceColor = [0.7 0.7 0.7];
plot([imuTimeS(start) imuTimeS(finish)],[thresholdTrans thresholdTrans], 'r:','LineWidth',2.4)

legend('DC-5hz','5-30hz','30-75hz','75-150hz','150-800hz','800-1200hz','1200> hz','threshold')
grid on
xlim([imuTimeS(start) imuTimeS(finish)])
% plot(imuTimeS(start:finish),totalVelocityTrans(start:finish),'g','LineWidth',3)

dpsViolationRate = size(find(totalVelocityTrans >= thresholdTrans),2)/L * 100
dpsViolationRate = size(find(totalVelocityTrans(start:finish) >= thresholdTrans),2)/(finish-start) * 100

count = 0;
for i= 1:L-40
if(~isempty(find(totalVelocityTrans(i:i+40) > thresholdTrans)))
count = count + 1;
end
end
count/(L-40) * 100

%% Total component plot

vel = [];
x_comp = [];
y_comp = [];
pixelPosAngle = pi/4;
pixelDist = 3;
for i = 1:length(freqL)
    x_comp(i,:) = filteredFreq_r(i,:) * 9 - pixelDist * filteredFreq_p(i,:) * cos(pixelPosAngle) + filteredFreq_v(i,:);
    y_comp(i,:) = filteredFreq_q(i,:) * 9 + pixelDist * filteredFreq_p(i,:) * sin(pixelPosAngle) - filteredFreq_w(i,:);
    vel(i,:) = sqrt(x_comp(i,:).^2 + y_comp(i,:).^2);
end

totalVelocity = sqrt((sum(x_comp)+0).^2 + (sum(y_comp)+0).^2);
threshold = 3.9; % m/s
residue = threshold - totalVelocity;
velBar = vertcat(vel,residue);


figure(11)
clf
subplot(2,1,1)
hold on
grid on
plot(imuTimeS(start:finish),totalVelocity(start:finish),'k','LineWidth',1)
plot([imuTimeS(start) imuTimeS(finish)],[threshold threshold], 'r--','LineWidth',2)
xlim([imuTimeS(start) imuTimeS(finish)])
title('Pointing error velocity - total');
xlabel('time[s]');
ylabel('m/s');
ylim([0 5])
legend('pointing velocity','threshold');
xlim([imuTimeS(start) imuTimeS(finish)])
% xlim([308 309])

subplot(2,1,2)
a = area(imuTimeS(start:finish),velBar(1:end-1,start:finish)');
xlim([imuTimeS(start) imuTimeS(finish)])
xlabel('time[s]')
ylabel('m/s')
title('Contribution ratio area plot - total')
hold on
plot([imuTimeS(start) imuTimeS(finish)],[threshold threshold], 'r:','LineWidth',2.4)
legend('DC-5hz','5-30hz','30-75hz','75-150hz','150-800hz','800-1200hz','1200> hz','threshold')
grid on
plot([],[],'m:','LineWidth',1.5)
plot([],[],'g:','LineWidth',1.5)
ylim([0 5])
xlim([imuTimeS(start) imuTimeS(finish)])
% xlim([308 309])

componentRatio = [];
meanVelBar = mean(velBar(1:end-1,start:finish),2);
stdVelBar = std(velBar(1:end-1,start:finish)');
totalVelBar = sum(meanVelBar);
componentRatio = meanVelBar/totalVelBar*100;

% figure(11111111)
% bar(componentRatio)
% xticks([1,2,3,4,5,6,7]);
% xticklabels({'1-5hz','5-30hz','30-75hz','75-150hz','150-800hz','800-1200hz','1200~hz'});
% ylabel('Contribution [%]')
% grid on
% 
% figure(11111112)
% bar(meanVelBar)
% xticks([1,2,3,4,5,6,7]);
% xticklabels({'1-5hz','5-30hz','30-75hz','75-150hz','150-800hz','800-1200hz','1200~hz'});
% ylabel('mean pointing velocity [m/s]')
% grid on


syncLen = 20;
count = 0;
missionLength = finish - start - syncLen;
for i = start:finish-syncLen
    if (~isempty(find(totalVelocity(i:i+syncLen)>threshold)))
        count = count + 1;
    end
end
totalSatis = count / (missionLength-syncLen) * 100;

totalSatis
meanTotVel = mean(totalVelocity(start:finish))
stdTotVel = std(totalVelocity(start:finish))

%% Contribution analysis
% save('oksang_wind_3','meanVelBar','stdVelBar','meanTotVel','stdTotVel')


