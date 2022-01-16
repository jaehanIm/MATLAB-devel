addpath('./..')

%% Load data

% 220112 inhouse - 10hz
% imuData = readtable('/home/jaehan/log/220112_[0]_inhouse/aSensorLog_220112_111012.csv');
% gdLog = readtable('/home/jaehan/log/220112_[0]_inhouse/gdLog_220112_111012.csv');

% 220112 inhouse - 20hz
% imuData = readtable('/home/jaehan/log/220112_[0]_inhouse/aSensorLog_220112_110312.csv');
% gdLog = readtable('/home/jaehan/log/220112_[0]_inhouse/gdLog_220112_110312.csv');

% 220113 oksang - 20hz
% imuData = readtable('/home/jaehan/log/220113_[0]_oksang/220113_090300/aSensorLog_220113_090300.csv');
% gdLog = readtable('/home/jaehan/log/220113_[0]_oksang/220113_090300/gdLog_220113_090300.csv');

% 220113 oksang - 10hz (noon)
% imuData = readtable('/home/jaehan/log/220113_[0]_oksang/220113_140901/aSensorLog_220113_140901.csv');
% gdLog = readtable('/home/jaehan/log/220113_[0]_oksang/220113_140901/gdLog_220113_140901.csv');

% 220113 oksang - 20hz (noon)
% imuData = readtable('/home/jaehan/log/220113_[0]_oksang/220113_141431/aSensorLog_220113_141431.csv');
% gdLog = readtable('/home/jaehan/log/220113_[0]_oksang/220113_141431/gdLog_220113_141431.csv');

% 220114 oksang - 10hz
imuData = readtable('/home/jaehan/log/220114_[0]_oksang/220114_132544/aSensorLog_220114_132544.csv');
gdLog = readtable('/home/jaehan/log/220114_[0]_oksang/220114_132544/gdLog_220114_132544.csv');

% 220114 oksang - 20hz
% imuData = readtable('/home/jaehan/log/220114_[0]_oksang/220114_133507/aSensorLog_220114_133507.csv');
% gdLog = readtable('/home/jaehan/log/220114_[0]_oksang/220114_133507/gdLog_220114_133507.csv');

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

% time synchronization
imuTime = datetime(imuData.rosTime,'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
gdTime = datetime(gdLog.rosTime(1:end),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
imuTimeDelay = seconds(imuTime(1) - gdTime(1));
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);
imuTimeS = seconds(imuTime - imuTime(1)) + imuTimeDelay;

% Assign data
acc_0 = imuData.acc_mpss_0;
acc_1 = imuData.acc_mpss_1;
acc_2 = imuData.acc_mpss_2;
gyro_0 = imuData.gyro_dps_0;
gyro_1 = imuData.gyro_dps_1;
gyro_2 = imuData.gyro_dps_2;

Fs = round(1/mean(diff(imuTimeS)));

disp("Data Loading Complete!")

%% pre processing
acc_0 = detrend(acc_0,1);
acc_1 = detrend(acc_1,1);
acc_2 = detrend(acc_2,2);
gyro_0 = detrend(gyro_0,1);
gyro_1 = detrend(gyro_1,1);
gyro_2 = detrend(gyro_2,1);

%% Raw plot
figure(1)
clf

subplot(6,1,1)
plot(imuTimeS,gyro_0)
title('gyro R')
xlabel('time[s]');
ylabel('[rad/s]');
hold on
ylim([-10*std(gyro_0) 10*std(gyro_0)])

subplot(6,1,2)
plot(imuTimeS,gyro_1)
title('gyro P')
xlabel('time[s]')
ylabel('[rad/s]');
hold on
ylim([-10*std(gyro_1) 10*std(gyro_1)])

subplot(6,1,3)
plot(imuTimeS,gyro_2)
title('gyro Y')
xlabel('time[s]')
ylabel('[rad/s]');
hold on
ylim([-10*std(gyro_2) 10*std(gyro_2)])

subplot(6,1,4)
plot(imuTimeS,acc_0)
title('acc X')
xlabel('time[s]');
ylabel('[m/s^2]');
hold on
ylim([-10*std(acc_0) 10*std(acc_0)])

subplot(6,1,5)
plot(imuTimeS,acc_1)
title('acc Y')
xlabel('time[s]')
ylabel('[m/s^2]');
hold on
ylim([-10*std(acc_1) 10*std(acc_1)])

subplot(6,1,6)
plot(imuTimeS,acc_2)
title('acc Z')
xlabel('time[s]')
ylabel('[m/s^2]');
hold on
ylim([-10*std(acc_2) 10*std(acc_2)])

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
plot(gdTimeS,gdLog.velUVW_mps_1)
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
    gduF = gdLog.velUVW_mps_0(i);
    gduB = gdLog.velUVW_mps_0(i+1);
    gdvF = gdLog.velUVW_mps_1(i);
    gdvB = gdLog.velUVW_mps_1(i+1);
    gdwF = gdLog.velUVW_mps_2(i);
    gdwB = gdLog.velUVW_mps_2(i+1);
    
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
plot(gdTimeS,gdLog.velUVW_mps_1)
title('velocity trend comparison')
legend('integrated v (HP filtered)','gdLog v')

disp('Partial FDI complete!')
%% data selection
y = gyro_0;

%% FFT plot
figure(2)
title('FFT')
clf
[freq, fftResult, psdResult] = data2fftpsd(u,Fs);
% plot(freq,fftResult);
plot(freq,fftResult);
hold on
[freq, fftResult, psdResult] = data2fftpsd(v,Fs);
% plot(freq,fftResult);
plot(freq,fftResult);
grid on
xlabel('hz');
ylabel('magnitude')
legend('u','v')
xlim([0.1 2000])

figure(3)
title('PSD')
clf
[freq, fftResult, psdResult] = data2fftpsd(gyro_0,Fs);
plot(freq,fftResult);
hold on
[freq, fftResult, psdResult] = data2fftpsd(gyro_1,Fs);
plot(freq,fftResult);
grid on
xlabel('hz');
ylabel('psd')
legend('gyro_R','gyro_P')
xlim([0.1 2000])

%% Spectrogram
y = gyro_0;
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

figure(8)
imagesc([freq(1) freq(end)],[imuTimeS(1) imuTimeS(end)],spectrogram(1:end,1:end));
xlabel('freq [hz]')
ylabel('time [s]')
% xlim([0 200])
% ylim([0 250])
title('Spectrogram Analysis [dt = 1.0s]')
colorbar

figure(88)
clf
sgtitle('spectrogram')

subplot(1,4,1)
imagesc([0 200],[imuTimeS(1) imuTimeS(end)],spectrogram(:,0*timeStep+1:200*timeStep));
colorbar
xlabel('hz')
ylabel('time [s]')
title('0~200hz')
hold on


subplot(1,4,2)
imagesc([40 500],[imuTimeS(1) imuTimeS(end)],spectrogram(:,40*timeStep:500*timeStep));
hold on
xlabel('hz')
ylabel('time [s]')
title('40~500hz')
colorbar

subplot(1,4,3)
imagesc([500 900],[imuTimeS(1) imuTimeS(end)],spectrogram(:,500*timeStep:900*timeStep));
hold on
xlabel('hz')
ylabel('time [s]')
title('500~900')
colorbar

subplot(1,4,4)
imagesc([900 1300],[imuTimeS(1) imuTimeS(end)],spectrogram(:,900*timeStep:1300*timeStep));
hold on
xlabel('hz')
ylabel('time [s]')
title('900~1300hz')
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
plot(gdTimeS,gdLog.pqr_dps_0)
plot(gdTimeS,gdLog.pqr_dps_1)
plot(gdTimeS,gdLog.pqr_dps_2)
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

freqL = [0 5  30  100 500 800 1200];
freqH = [5 30 100 500 800 1200 2000];
% freqL = [0 5  30  100 500 1200];
% freqH = [5 30 100 500 800 2000];
% freqL = [0 5  15  100 500 1200];
% freqH = [5 15 25 500 800 2000];
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

freqL = [0 5  30  100 500 1200];
freqH = [5 30 100 500 800 2000];
% freqL = [0 5  30  100 300 500 1200];
% freqH = [5 30 100 300 500 800 2000];
% freqL = [0 5  15  100 500 1200];
% freqH = [5 15 25 500 800 2000];
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

%% vibration component plot

start = 1;
finish = L;
% start = imuJobStartIdx(2);
% finish = imuJobEndIdx(2);

figure(9)
clf
hold on
grid on
plot(imuTimeS(start:finish),totalVelocity(start:finish),'k','LineWidth',1)
% for i = 1:size(vel,1)
% plot(imuTimeS(start:finish),vel(i,start:finish))
% end
plot([imuTimeS(start) imuTimeS(finish)],[24.8 24.8], 'r--','LineWidth',2)
xlim([imuTimeS(start) imuTimeS(finish)])
title('Pointing error velocity - angular component');
legend('pointing velocity','threshold');
xlabel('time[s]');
ylabel('deg/s');
% xlim([110 110.1])


figure(99)
clf
a = area(imuTimeS(start:finish),velBar(1:end-1,start:finish)');
xlim([imuTimeS(start) imuTimeS(finish)])
xlabel('time[s]')
ylabel('deg/s')
title('Contribution ratio area plot - angular component')
hold on
% a(8).FaceColor = [0.7 0.7 0.7];
plot([imuTimeS(start) imuTimeS(finish)],[24.8 24.8], 'r:','LineWidth',2.4)
legend('DC-5hz','5-30hz','30-100hz','100-500hz','500-800hz','800-1200hz','1200~hz','residue','threshold')
grid on
% xlim([110 110.1])
% plot(imuTimeS(start:finish),totalVelocity(start:finish),'g','LineWidth',3)

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

% start = imuJobStartIdx(2);
% finish = imuJobEndIdx(2);

figure(10)
clf
hold on
grid on
plot(imuTimeS(start:finish),totalVelocityTrans(start:finish),'k','LineWidth',1)
plot([imuTimeS(start) imuTimeS(finish)],[thresholdTrans thresholdTrans], 'r--','LineWidth',2)
xlim([imuTimeS(start) imuTimeS(finish)])
title('Pointing error velocity - translational component');
legend('pointing velocity','threshold');
xlabel('time[s]');
ylabel('m/s');
% xlim([110 110.1])
figure(110)
clf
a = area(imuTimeS(start:finish),velBarTrans(1:end-1,start:finish)');
xlim([imuTimeS(start) imuTimeS(finish)])
xlabel('time[s]')
ylabel('m/s')
title('Contribution ratio area plot - translational component')
hold on
% a(8).FaceColor = [0.7 0.7 0.7];
plot([imuTimeS(start) imuTimeS(finish)],[thresholdTrans thresholdTrans], 'r:','LineWidth',2.4)

legend('1-5hz','5-30hz','30-100hz','100-500hz','500-800hz','800-1200hz','1200~hz','residue','threshold')
grid on
% xlim([110 110.1])
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

start = 600000;
finish = 600400;

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

start = 1;
finish = L;

figure(11)
clf
hold on
grid on
plot(imuTimeS(start:finish),totalVelocity(start:finish),'k','LineWidth',1)
plot([imuTimeS(start) imuTimeS(finish)],[threshold threshold], 'r--','LineWidth',2)
xlim([imuTimeS(start) imuTimeS(finish)])
title('Pointing error velocity - total');
xlabel('time[s]');
ylabel('m/s');
% xlim([110 110.5])

legend('pointing velocity','threshold');
ylim([0 40])

figure(1111)
clf
hold on
grid on
plot(totalVelocity(start:finish),'k','LineWidth',1)
title('Pointing error velocity - total'); 
xlabel('time[s]');
ylabel('m/s');
% xlim([110 110.5])

legend('pointing velocity','threshold');
ylim([0 40])

%% Total contribution

move = -400 * 0;
start = 560000 + move;
finish = 560000 + move + 200;

figure(111111)
clf
a = area(imuTimeS(start:finish),velBar(1:end-1,start:finish)');
xlim([imuTimeS(start) imuTimeS(finish)])
xlabel('time[s]')
ylabel('m/s')
title('Contribution ratio area plot - total')
hold on
% a(8).FaceColor = [0.7 0.7 0.7];
plot([imuTimeS(start) imuTimeS(finish)],[threshold threshold], 'r:','LineWidth',2.4)
legend('1-5hz','5-30hz','30-100hz','100-500hz','500-800hz','800-1200hz','1200~hz','residue','threshold')
grid on
% xlim([110 110.5])
% plot(imuTimeS(start:finish),totalVelocity(start:finish),'g','LineWidth',3)
plot([],[],'m:','LineWidth',1.5)
plot([],[],'g:','LineWidth',1.5)


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
mean(totalVelocity(start:finish))

%% Contribution analysis
% temp = sum(velBar(:,6000:10000),2)/4000
% save('10hz_inhouse.mat','temp')

% temp = sum(velBar(:,460000:464000),2)/4000
% save('0hz_inhouse.mat','temp')

% temp = sum(velBar(:,6000:10000),2)/4000
% save('20hz_inhouse.mat','temp')

% temp = sum(velBar(:,start:start+4000),2)/4000;
% save('10hz_210113_1.mat','temp');

% temp = sum(velBar(:,560000:564000),2)/4000;
% save('10hz_220114.mat','temp');
% 
% temp = sum(velBar(:,680000:684000),2)/4000;
% save('20hz_220114.mat','temp');