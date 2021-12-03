% Load data
% imuData = readtable('/home/jaehan/log/rooftop_vibration/211124_160250_rooftop_20Hz_noPic/aSensorLog_211124_160250.csv');
% gdLog = readtable('/home/jaehan/log/rooftop_vibration/211124_160250_rooftop_20Hz_noPic/gdLog_211124_160250.csv');
% imuData = readtable('/home/jaehan/log/rooftop_vibration/211124_162252_rooftop_10Hz_Pic/aSensorLog_211124_162252.csv');
% gdLog = readtable('/home/jaehan/log/rooftop_vibration/211124_162252_rooftop_10Hz_Pic/gdLog_211124_162252.csv');
% imuData = readtable('/home/jaehan/log/rooftop_vibration/211124_161133_rooftop_20Hz_Pic/aSensorLog_211124_161133.csv');
% gdLog = readtable('/home/jaehan/log/rooftop_vibration/211124_161133_rooftop_20Hz_Pic/gdLog_211124_161133.csv');


jobIdx = [457873 707391 883091 1053720 1181680 1470740 1646400];
% last two = start and end of hold test
parseStart = jobIdx(1);
parseEnd = jobIdx(7);

% parseStart = 800000;
% parseEnd = 1400000;
% 
% parseStart = 816600; %302419;
% parseEnd = 830600;%1743950;
% parseStart = 302419;
% parseEnd = 1743950;
% parseStart = 300000;
% parseEnd = 340000;

imuTimeS = imuData.rosTime(parseStart:parseEnd) - imuData.rosTime(1);
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);
imuTime = datetime(imuData.rosTime(parseStart:parseEnd),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
gdTime = datetime(gdLog.rosTime(1:end),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');


acc_0 = imuData.acc_mpss_0(parseStart:parseEnd);
acc_1 = imuData.acc_mpss_1(parseStart:parseEnd);
acc_2 = imuData.acc_mpss_2(parseStart:parseEnd);
gyro_0 = imuData.gyro_dps_0(parseStart:parseEnd);
gyro_1 = imuData.gyro_dps_1(parseStart:parseEnd);
gyro_2 = imuData.gyro_dps_2(parseStart:parseEnd);

Fs = round(1/mean(diff(imuTimeS)));

%% pre processing
acc_0 = detrend(acc_0,1);
acc_1 = detrend(acc_1,1);
acc_2 = detrend(acc_2,2);
gyro_0 = detrend(gyro_0,1);
gyro_1 = detrend(gyro_1,1);
gyro_2 = detrend(gyro_2,1);

%% Window analysis

%% Raw plot
figure(1)
clf

subplot(3,1,1)
plot(imuTimeS,gyro_0)
title('gyro R')
xlabel('time[s]');
ylabel('[rad/s]');
hold on
plot([(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],[-0.6 0.6],'k--')

subplot(3,1,2)
plot(imuTimeS,gyro_1)
title('gyro P')
xlabel('time[s]')
ylabel('[rad/s]');
hold on
plot([(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],[-0.5 0.5],'k--')

subplot(3,1,3)
plot(imuTimeS,gyro_2)
title('gyro Y')
xlabel('time[s]')
ylabel('[rad/s]');
hold on
plot([(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],[-2 2],'k--')

figure(11)
clf
subplot(3,1,1)
plot(imuTimeS,acc_0)
title('acc X')
xlabel('time[s]');
ylabel('[m/s^2]');
hold on
plot([(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],[-0.6 0.6],'k--')
plot([(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],[-0.6 0.6],'k--')

subplot(3,1,2)
plot(imuTimeS,acc_1)
title('acc Y')
xlabel('time[s]')
ylabel('[m/s^2]');
hold on
plot([(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],[-0.5 0.5],'k--')
plot([(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],[-0.5 0.5],'k--')

subplot(3,1,3)
plot(imuTimeS,acc_2)
title('acc Z')
xlabel('time[s]')
ylabel('[m/s^2]');
hold on
plot([(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],[-2 2],'k--')
plot([(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],[-2 2],'k--')

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
filterFreqUpper = 0.15; %hz
U = fft(u);
V = fft(v);
W = fft(w);

df = L/Fs;
kernel = ones(L,1);
kernel(1:filterFreqUpper*df) = 0;
kernel(L-filterFreqUpper*df:end) = 0;

U_filter = U .* kernel;
V_filter = V .* kernel;
W_filter = W .* kernel;

u = real(ifft(U_filter));
v = real(ifft(V_filter));
w = real(ifft(W_filter));

figure(4)
clf
hold on
plot(imuTimeS,w)
% plot(imuTimeS,filteredv)
% plot(imuTimeS,detrend(filteredv,0))
grid on
xlabel('time[s]')
ylabel('v [m/s]');
plot(gdTimeS-80,gdLog.velUVW_mps_2)
legend('integrated v (HP filtered)','gdLog v')
title('velocity trend comparison')
xlim([10 60])
ylim([-5 5])


%% data selection
y = v;

%% data selection 2
% y = v' + v(1);

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

%% FFT decomposition
figure(7)
clf
hold on
grid on
Y = fft(y);

freqL = [1 5  30  100 500 800 1200];
freqH = [5 30 100 500 800 1200 2000];
df = L/Fs;
filteredFreq_R = [];

for i = 1:length(freqL)
    P_filter = Y;
    kernel = zeros(L,1);
    kernel(freqL(i)*df:freqH(i)*df) = 1;
    kernel(L-freqH(i)*df:L-freqL(i)*df) = 1;
    P_filter = P_filter .* kernel;
    fftRegen = ifft(P_filter);
    
    filteredFreq_R(i,:) = real(fftRegen);
    
    subplot(length(freqL),1,i)
    plot(imuTimeS,real(fftRegen))
    xlim([10 11])
    labelText = num2str(freqL(i)) + "hz to " + num2str(freqH(i)) + "hz";
    title(labelText);
end

sgtitle('FFT spectral decompostion','fontsize',14) 

%% EMD decomposition
% IMFresult = EMD(imuTimeS, y);
% figure(5)
% clf
% hold on
% grid on
% plot(imuTimeS,y)
% for i = 1:8
%     subplot(9,1,i+1)
%     plot(imuTimeS, IMFresult(i,:))
%     xlim([11 12])
%     titleText = "IMF " + num2str(i);
%     title(titleText)
% end
% sgtitle('IMF decomposition') 

%% SVD decomposition

ySel = y(1:Fs);
LSel = length(ySel);
colNum = 100;
hankel = zeros(LSel+1-colNum,colNum);
for i = 1:LSel + 1 - colNum
    hankel(i,:) = ySel(i:i+colNum-1);
end

[U,S,V] = svd(hankel);
test_regen_full = U * S * V';
% rate = 10;
% test_regen = U(:,1:rate) * S(1:rate,1:rate) * V(:,1:rate)';
% 
% signal_regen_r = test_regen(1,:);
% signal_regen_c = test_regen(2:end,end)';
% signal_regen = [signal_regen_r,signal_regen_c];

svdResult = [];

figure(6)
clf
hold on
grid on
for i = 1:20
    test_regen = U(:,i) * S(i,i) * V(:,i)';
    signal_regen_r = test_regen(1,:);
    signal_regen_c = test_regen(2:end,end)';
    signal_regen = [signal_regen_r,signal_regen_c];
    subplot(20,1,i)
    svdResult(i,:) = signal_regen;
    plot(imuTimeS(1:Fs),signal_regen)
%     ylim([-0.1 0.1])
end

sgtitle('SVD decomposition [3s duration]') 

%% Spectrogram

timeStep = .1; %[s]

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
imagesc([freq(1) freq(end)],[imuTimeS(1) imuTimeS(end)],spectrogram(:,1:end));
xlabel('freq [hz]')
ylabel('time [s]')
% xlim([0 200])
% ylim([0 250])
title('Spectrogram Analysis [dt = 1.0s]')
colorbar

figure(88)
clf
sgtitle('spectrogram')

subplot(1,3,1)
imagesc([0 200],[imuTimeS(1) imuTimeS(end)],spectrogram(:,1:200*timeStep));
colorbar
xlabel('hz')
ylabel('time [s]')
title('0~200hz')
hold on
plot([0 200],[(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([0 200],[(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([0 200],[(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([0 200],[(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([0 200],[(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([0 200],[(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],'m--','LineWidth',1.5)

subplot(1,3,2)
imagesc([100 500],[imuTimeS(1) imuTimeS(end)],spectrogram(:,100*timeStep:500*timeStep));
hold on
plot([100 500],[(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([100 500],[(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([100 500],[(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([100 500],[(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([100 500],[(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([100 500],[(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],'m--','LineWidth',1.5)
xlabel('hz')
ylabel('time [s]')
title('100~500hz')
colorbar

subplot(1,3,3)
imagesc([900 1300],[imuTimeS(1) imuTimeS(end)],spectrogram(:,900*timeStep:1300*timeStep));
hold on
plot([900 1300],[(jobIdx(1)-parseStart)/Fs (jobIdx(1)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([900 1300],[(jobIdx(2)-parseStart)/Fs (jobIdx(2)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([900 1300],[(jobIdx(3)-parseStart)/Fs (jobIdx(3)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([900 1300],[(jobIdx(4)-parseStart)/Fs (jobIdx(4)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([900 1300],[(jobIdx(5)-parseStart)/Fs (jobIdx(5)-parseStart)/Fs],'m--','LineWidth',1.5)
plot([900 1300],[(jobIdx(6)-parseStart)/Fs (jobIdx(6)-parseStart)/Fs],'m--','LineWidth',1.5)
xlabel('hz')
ylabel('time [s]')
title('900~1300hz')
colorbar

%% Vibration component analysis - pqr
p = gyro_0;
q = gyro_1;
r = gyro_2;

P = fft(p);
Q = fft(q);
R = fft(r);

freqL = [1 5  30  100 500 800 1200];
freqH = [5 30 100 500 800 1200 2000];
df = L/Fs;

filteredFreq_p = [];
filteredFreq_q = [];
filteredFreq_r = [];

for i = 1:length(freqL)
    P_filter = P;
    Q_filter = Q;
    R_filter = R;
    
    kernel = zeros(L,1);
    kernel(freqL(i)*df:freqH(i)*df) = 1;
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

freqL = [1 5  30  100 500 800 1200];
freqH = [5 30 100 500 800 1200 2000];
df = L/Fs;

filteredFreq_u = [];
filteredFreq_v = [];
filteredFreq_w = [];

for i = 1:length(freqL)
    U_filter = U;
    V_filter = V;
    W_filter = W;
    
    kernel = zeros(L,1);
    kernel(freqL(i)*df:freqH(i)*df) = 1;
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
legend('1-5hz','5-30hz','30-100hz','100-500hz','500-800hz','800-1200hz','1200~hz','residue','threshold')
grid on
% xlim([110 110.1])
% plot(imuTimeS(start:finish),totalVelocity(start:finish),'g','LineWidth',3)

dpsViolationRate = size(find(totalVelocity >= threshold),2)/L * 100
dpsViolationRate = size(find(totalVelocity(start:finish) >= threshold),2)/(finish-start) * 100

count = 0;
for i= 1:L-40
if(~isempty(find(totalVelocity(i:i+40) > threshold)))
count = count + 1;
end
end
count/(L-40) * 100

%% vibration component plot translational

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
vel = [];
x_comp = [];
y_comp = [];
pixelPosAngle = pi/2;
pixelDist = 0;
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
hold on
grid on
plot(imuTimeS(start:finish),totalVelocity(start:finish),'k','LineWidth',1)
plot([imuTimeS(start) imuTimeS(finish)],[threshold threshold], 'r--','LineWidth',2)
xlim([imuTimeS(start) imuTimeS(finish)])
title('Pointing error velocity - total');
legend('pointing velocity','threshold');
xlabel('time[s]');
ylabel('m/s');
% xlim([110 110.5])

figure(1111)
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


dpsViolationRate = size(find(totalVelocity(start:finish) >= threshold),2)/(finish-start) * 100

count = 0;
for i= 1:L-40
if(~isempty(find(totalVelocity(i:i+40) > threshold)))
count = count + 1;
end
end
count/(L-40) * 100