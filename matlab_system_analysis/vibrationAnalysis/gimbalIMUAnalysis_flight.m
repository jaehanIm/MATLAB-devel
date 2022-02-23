addpath('./..')
%% Load data

% imuData = readtable('/home/jaehan/log/211126/211126_110526/aSensorLog_211126_110526.csv');
% gdLog = readtable('/home/jaehan/log/211126/211126_110526/gdLog_211126_110526.csv');
% photoDir = dir('/home/jaehan/Desktop/M300G_Nikon_Test/1');
% photoXml = readtable('/home/jaehan/Desktop/M300G_Nikon_Test/df_xmp_1.csv');
% faultPhotoIdx = [10 13 15 17 29 43 65 67 68 69 70 71 72 73 74 75 79 81 82 85 86 87 91 97 104 110 125 126];
% realBadPhotoIdx = [15    29    65    67    70    74    81    85    86   110];
% photoTimeTrim = seconds(0.8);

% imuData = readtable('/home/jaehan/log/211126/211126_113400/aSensorLog_211126_113400.csv');
% gdLog = readtable('/home/jaehan/log/211126/211126_113400/gdLog_211126_113400.csv');
% photoDir = dir('/home/jaehan/Desktop/M300G_Nikon_Test/2');
% photoXml = readtable('/home/jaehan/Desktop/M300G_Nikon_Test/df_xmp_2.csv');
% faultPhotoIdx = [26 44 69 85 88 93 113 128 129];
% realBadPhotoIdx = [69];
% photoTimeTrim = seconds(1.0);

% 220210 ss1
% imuData = readtable('/home/jaehan/log/220210_temp/220210_141454_soso1/aSensorImu_220210_141454.csv');
% imuData(1:165,:) = [];
% % imuData(1680000:end,:) = []; % 6-1
% imuData(1:1680000,:) = []; % 6-2
% gdLog = readtable('/home/jaehan/log/220210_temp/220210_141454_soso1/gdLog_220210_141454.csv');
% % gdLog(21600:end,:) = []; % 6-1
% gdLog(1:21600,:) = []; % 6-2
% photoDir = dir('/home/jaehan/log/220210_temp/220210_141454_soso1/soft-soft sortie 1');
% photoXml = readtable('/home/jaehan/log/220210_temp/220210_141454_soso1/df_xmp.csv');
% photoTimeTrim = seconds(0.8 + 0.6780 - 0.1);

% 220210 ss2
% imuData = readtable('/home/jaehan/log/220210_temp/220210_143844_soso2/aSensorImu_220210_143844.csv');
% imuData(1:165,:) = [];
% % imuData(1:1978220-160000,:) = []; % 6-1
% imuData(1978220:end,:) = [];  % 6-2
% gdLog = readtable('/home/jaehan/log/220210_temp/220210_143844_soso2/gdLog_220210_143844.csv');
% % gdLog(1:22542,:) = [];
% gdLog(22542:end,:) = [];
% photoDir = dir('/home/jaehan/log/220210_temp/220210_143844_soso2/soft-soft sortie 2');
% photoXml = readtable('/home/jaehan/log/220210_temp/220210_143844_soso2/df_xmp.csv');
% photoTimeTrim = seconds(0.8 + 0.6780 - 0.1 + 0.0730);

% 220210 st1
% imuData = readtable('/home/jaehan/log/220210_temp/220210_150409_stso1/aSensorImu_220210_150409.csv');
% imuData(1:165,:) = [];
% % imuData(1731564:end,:) = []; % 6-1
% imuData(1:1731564,:) = []; % 6-2
% gdLog = readtable('/home/jaehan/log/220210_temp/220210_150409_stso1/gdLog_220210_150409.csv');
% % gdLog(19179:end,:) = [];
% gdLog(1:19179,:) = [];
% photoDir = dir('/home/jaehan/log/220210_temp/220210_150409_stso1/stiff-soft sortie 1');
% photoXml = readtable('/home/jaehan/log/220210_temp/220210_150409_stso1/df_xmp.csv');
% photoTimeTrim = seconds(0.8 + 0.6780 - 0.1);

% 220210 st2
% imuData = readtable('/home/jaehan/log/220210_temp/220210_163334_stso2/aSensorImu_220210_163334.csv');
% imuData(1:165,:) = [];
% % imuData(1:1731564,:) = [];
% imuData(1731564:end,:) = [];
% gdLog = readtable('/home/jaehan/log/220210_temp/220210_163334_stso2/gdLog_220210_163334.csv');
% % gdLog(1:20559,:) = [];
% gdLog(20559:end,:) = [];
% photoDir = dir('/home/jaehan/log/220210_temp/220210_163334_stso2/stiff-soft sortie 4');
% photoXml = readtable('/home/jaehan/log/220210_temp/220210_163334_stso2/df_xmp.csv');
% photoTimeTrim = seconds(0.8 + 0.6780 - 0.1550);

% 220211 ss1



%% Preprocessing

photoDir(1:2) = [];
photoDir(end) = [];
photoL = length(photoDir); 

photoTime = photoXml.xmpTime;
photoTime = datetime(photoTime,'TimeZone','Asia/Tokyo');
photoTime = photoTime + photoTimeTrim;


% Job index parser
jobStartIdx = [];
jobEndIdx = [];
landIdx = [];
for i = 2:size(gdLog,1)
    if gdLog.fcMcMode(i) == 2 && gdLog.fcMcMode(i-1) == 1
        jobStartIdx = horzcat(jobStartIdx,i);
    elseif gdLog.fcMcMode(i) == 1 && gdLog.fcMcMode(i-1) == 2
        jobEndIdx = horzcat(jobEndIdx,i);
    elseif gdLog.fcMcMode(i) == 255 && gdLog.fcMcMode(i-1) == 0
        landIdx = i;
    end
end


% jobName = ["SME","#2 LE bf","PP","#2 UPSS bf","PP","#1 LTSS bf","PP","#1 TE bf","PP","#0 UPSS bf","PP","#0 TE bf", "ShutOff"];
jobName = ["SME","PP","#0 DNSS","PP","#0 LE","PP","#1 RTSS","PP","#1 LE","PP","#2 DNSS","PP","#2 TE","ShutOff"];
% jobName6_2 = ["SME","","","",""];
parseStart = 1;
parseEnd = size(imuData,1);

% time synchronization
imuTime = datetime(imuData.rosTime(parseStart:parseEnd),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
gdTime = datetime(gdLog.rosTime(1:end),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
imuTimeDelay = seconds(imuTime(1) - gdTime(1));
photoTimeDelay = photoTime(1) - gdTime(1);
imuTimeS = imuData.rosTime(parseStart:parseEnd) - imuData.rosTime(1) + imuTimeDelay;
photoTimeS = seconds(photoTime - photoTime(1) + photoTimeDelay);
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);
photoTimeS(photoTimeS<=0) = 0;

landIdx = 10000;
jobStartIdxS = gdTimeS(jobStartIdx);
jobEndIdxS = gdTimeS(jobEndIdx);
landIdxS = gdTimeS(landIdx);

% Job index to imudata idx
imuJobStartIdx = [];
imuJobEndIdx = [];
for i = 1:length(jobStartIdx)
    imuJobStartIdx(i) = findNearestIdx(imuTimeS,gdTimeS(jobStartIdx(i)));
    imuJobEndIdx(i) = findNearestIdx(imuTimeS,gdTimeS(jobEndIdx(i)));
end

photoJobCategory = [];
for i = 1:length(jobStartIdx)
    photoJobCategory{i} = find(photoTimeS > jobStartIdxS(i) & photoTimeS < jobEndIdxS(i));
end

% Assign data
acc_0 = imuData.acc_mpss_0(parseStart:parseEnd);
acc_1 = imuData.acc_mpss_1(parseStart:parseEnd);
acc_2 = imuData.acc_mpss_2(parseStart:parseEnd);
gyro_0 = imuData.gyro_dps_0(parseStart:parseEnd);
gyro_1 = imuData.gyro_dps_1(parseStart:parseEnd);
gyro_2 = imuData.gyro_dps_2(parseStart:parseEnd);

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
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(gyro_0)*5 std(gyro_0)*5],'r--')
    plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(gyro_0)*5 std(gyro_0)*5],'g--')
end
ylim([-10*std(gyro_0) 10*std(gyro_0)])

subplot(6,1,2)
plot(imuTimeS,gyro_1)
title('gyro P')
xlabel('time[s]')
ylabel('[rad/s]');
hold on
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(gyro_1)*5 std(gyro_1)*5],'r--')
    plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(gyro_1)*5 std(gyro_1)*5],'g--')
end
ylim([-10*std(gyro_1) 10*std(gyro_1)])

subplot(6,1,3)
plot(imuTimeS,gyro_2)
title('gyro Y')
xlabel('time[s]')
ylabel('[rad/s]');
hold on
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(gyro_2)*5 std(gyro_2)*5],'r--')
    plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(gyro_2)*5 std(gyro_2)*5],'g--')
end
ylim([-10*std(gyro_2) 10*std(gyro_2)])

subplot(6,1,4)
plot(imuTimeS,acc_0)
title('acc X')
xlabel('time[s]');
ylabel('[m/s^2]');
hold on
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(acc_0)*5 std(acc_0)*5],'r--')
    plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(acc_0)*5 std(acc_0)*5],'g--')
end
ylim([-10*std(acc_0) 10*std(acc_0)])

subplot(6,1,5)
plot(imuTimeS,acc_1)
title('acc Y')
xlabel('time[s]')
ylabel('[m/s^2]');
hold on
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(acc_1)*5 std(acc_1)*5],'r--')
    plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(acc_1)*5 std(acc_1)*5],'g--')
end
ylim([-10*std(acc_1) 10*std(acc_1)])

subplot(6,1,6)
plot(imuTimeS,acc_2)
title('acc Z')
xlabel('time[s]')
ylabel('[m/s^2]');
hold on
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(acc_2)*5 std(acc_2)*5],'r--')
    plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(acc_2)*5 std(acc_2)*5],'g--')
end
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
for i = 1:length(jobStartIdxS)
plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(v)*5 std(v)*5],'r--')
plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(v)*5 std(v)*5],'g--')
end
legend('integrated v (HP filtered)','gdLog v')


%% Partial FDI
for i = 1:length(gdTimeS)-1
    interest = find(imuTimeS>gdTimeS(i) & imuTimeS<gdTimeS(i+1));
    
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
for i = 1:length(jobStartIdxS)
plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(v)*5 std(v)*5],'r--')
plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(v)*5 std(v)*5],'g--')
end
legend('integrated v (HP filtered)','gdLog v')

disp('Partial FDI complete!')
%% data selection
y = gyro_0;

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
%     xlim([10 11])
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

% ySel = y(1:Fs);
% LSel = length(ySel);
% colNum = 100;
% hankel = zeros(LSel+1-colNum,colNum);
% for i = 1:LSel + 1 - colNum
%     hankel(i,:) = ySel(i:i+colNum-1);
% end
% 
% [U,S,V] = svd(hankel);
% test_regen_full = U * S * V';
% % rate = 10;
% % test_regen = U(:,1:rate) * S(1:rate,1:rate) * V(:,1:rate)';
% % 
% % signal_regen_r = test_regen(1,:);
% % signal_regen_c = test_regen(2:end,end)';
% % signal_regen = [signal_regen_r,signal_regen_c];
% 
% svdResult = [];
% 
% figure(6)
% clf
% hold on
% grid on
% for i = 1:20
%     test_regen = U(:,i) * S(i,i) * V(:,i)';
%     signal_regen_r = test_regen(1,:);
%     signal_regen_c = test_regen(2:end,end)';
%     signal_regen = [signal_regen_r,signal_regen_c];
%     subplot(20,1,i)
%     svdResult(i,:) = signal_regen;
%     plot(imuTimeS(1:Fs),signal_regen)
% %     ylim([-0.1 0.1])
% end
% 
% sgtitle('SVD decomposition [3s duration]') 

%% Spectrogram
y = gyro_0;
timeStep = .5; %[s]

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

subplot(1,4,1)
imagesc([0 200],[imuTimeS(1) imuTimeS(end)],spectrogram(:,0*timeStep+1:200*timeStep));
colorbar
xlabel('hz')
ylabel('time [s]')
title('0~200hz')
hold on
for i = 1:length(jobStartIdxS)
    plot([0 200],[jobStartIdxS(i) jobStartIdxS(i)],'m--','LineWidth',1.5)
    plot([0 200],[jobEndIdxS(i) jobEndIdxS(i)],'g--','LineWidth',1.5)
    text(2,jobStartIdxS(i)+2,jobName(i),'Color',[1 1 1])
end
plot([0 200], [landIdxS landIdxS], 'k--','LineWidth',1.5)
text(2,landIdxS+2,jobName(end),'Color',[1 1 1])


subplot(1,4,2)
imagesc([40 500],[imuTimeS(1) imuTimeS(end)],spectrogram(:,40*timeStep:500*timeStep));
hold on
for i = 1:length(jobStartIdxS)
    plot([40 500],[jobStartIdxS(i) jobStartIdxS(i)],'m--','LineWidth',1.5)
    plot([40 500],[jobEndIdxS(i) jobEndIdxS(i)],'g--','LineWidth',1.5)
end
plot([40 500], [landIdxS landIdxS], 'k--','LineWidth',1.5)
xlabel('hz')
ylabel('time [s]')
title('40~500hz')
colorbar

subplot(1,4,3)
imagesc([500 900],[imuTimeS(1) imuTimeS(end)],spectrogram(:,500*timeStep:900*timeStep));
hold on
for i = 1:length(jobStartIdxS)
    plot([500 900],[jobStartIdxS(i) jobStartIdxS(i)],'m--','LineWidth',1.5)
    plot([500 900],[jobEndIdxS(i) jobEndIdxS(i)],'g--','LineWidth',1.5)
end
plot([500 900], [landIdxS landIdxS], 'k--','LineWidth',1.5)
xlabel('hz')
ylabel('time [s]')
title('500~900')
colorbar

subplot(1,4,4)
imagesc([900 1300],[imuTimeS(1) imuTimeS(end)],spectrogram(:,900*timeStep:1300*timeStep));
hold on
for i = 1:length(jobStartIdxS)
    plot([900 1300],[jobStartIdxS(i) jobStartIdxS(i)],'m--','LineWidth',1.5)
    plot([900 1300],[jobEndIdxS(i) jobEndIdxS(i)],'g--','LineWidth',1.5)
end
plot([900 1300], [landIdxS landIdxS], 'k--','LineWidth',1.5)
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
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[0 200],'m--','LineWidth',1.5)
    plot([jobEndIdxS(i) jobEndIdxS(i)],[0 200],'g--','LineWidth',1.5)
    text(jobStartIdxS(i)+2,2,jobName(i),'Color',[1 1 1])
end
plot([landIdxS landIdxS], [0 200], 'k--','LineWidth',1.5)
text(landIdxS+2,2,jobName(end),'Color',[1 1 1])
xlim([imuTimeS(1) imuTimeS(end)])
ylim([0 200])

subplot(3,1,2)
hold on
grid on
plot(gdTimeS,gdLog.pqrBody_dps_0)
plot(gdTimeS,gdLog.pqrBody_dps_1)
plot(gdTimeS,gdLog.pqrBody_dps_2)
xlim([imuTimeS(1) imuTimeS(end)])
colorbar
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[-100 100],'m--','LineWidth',1.5)
    plot([jobEndIdxS(i) jobEndIdxS(i)],[-100 100],'g--','LineWidth',1.5)
    text(jobStartIdxS(i)+2,2,jobName(i),'Color',[1 1 1])
end
plot([landIdxS landIdxS], [-100 100], 'k--','LineWidth',1.5)
text(landIdxS+2,2,jobName(end),'Color',[1 1 1])

disp("Spectrogram plotting complete!")


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

freqL = [0.05 5  30  100 500 800 1200];
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
% start = imuJobStartIdx(4);
% finish = imuJobEndIdx(4);

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
for i= 1:L-20
if(~isempty(find(totalVelocity(i:i+20) > threshold)))
count = count + 1;
end
end
count/(L-20) * 100

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
for i = 1:length(jobStartIdxS)
plot([jobStartIdxS(i) jobStartIdxS(i)],[0 std(v)*5],'r--')
plot([jobEndIdxS(i) jobEndIdxS(i)],[0 std(v)*5],'g--')
end
for i = 1:length(jobStartIdx)
    text(jobStartIdxS(i)+2,3,jobName(i),'Fontsize',10,'FontWeight','bold','Color','m')
end
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
postShootCut = 0.34;
cutEnd = photoTimeS + postShootCut;
syncLen = 20;

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
hold on
grid on
plot(imuTimeS(start:finish),totalVelocity(start:finish),'k','LineWidth',1)
plot([imuTimeS(start) imuTimeS(finish)],[threshold threshold], 'r--','LineWidth',2)
xlim([imuTimeS(start) imuTimeS(finish)])
title('Pointing error velocity - total');
xlabel('time[s]');
ylabel('m/s');
% xlim([110 110.5])
for i = 1:size(jobStartIdx,2)
    plot([jobStartIdxS(i) jobStartIdxS(i)], [0 max(totalVelocity)],'r--')
    plot([jobEndIdxS(i) jobEndIdxS(i)], [0 max(totalVelocity)],'g--')
end
plot([landIdxS landIdxS], [0 max(totalVelocity)],'k--')

for i = 1:length(photoTimeS)
    plot(photoTimeS(i),0,'g*','LineWidth',4,'MarkerSize',10)
    plot(photoTimeS(i)+postShootCut,0,'r*','LineWidth',2,'MarkerSize',5)
end
% for i = 1:length(faultPhotoIdx)
%     plot(photoTimeS(faultPhotoIdx(i)),0,'rx','LineWidth',4,'MarkerSize',10)
% end
% for i = 1:length(realBadPhotoIdx)
%     plot(photoTimeS(realBadPhotoIdx(i)),0,'mo','LineWidth',4,'MarkerSize',10)
% end

for i = 1:length(jobStartIdx)
    text(jobStartIdxS(i)+2,30,jobName(i),'Fontsize',10,'FontWeight','bold','Color','m')
end

legend('pointing velocity','threshold');
ylim([0 40])


%% Total contribution

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
plot([],[],'m:','LineWidth',1.5)
plot([],[],'g:','LineWidth',1.5)


dpsViolationRate = size(find(totalVelocity(start:finish) >= threshold),2)/(finish-start) * 100

count = 0;
for i= 1:L-20
if(~isempty(find(totalVelocity(i:i+20) > threshold)))
count = count + 1;
end
end

shotNum = [];
for i = 1:length(jobStartIdx)
    shotNum = vertcat(shotNum,length(photoJobCategory{i}));
end

disp('violation rate')



%% Partial satisfaction rate

partialSatis = [];
for i = 1:length(imuJobStartIdx)
    jobLength = imuJobEndIdx(i) - imuJobStartIdx(i);
    count = 0;
    for j = imuJobStartIdx(i):imuJobEndIdx(i)
        if (~isempty(find(totalVelocity(j:j+20)>threshold)))
            if isempty(find((photoTimeS - imuTimeS(i)) .* (cutEnd - imuTimeS(i))<0))
                count = count + 1;
            end
        end
    end
    partialSatis(i) = count/(jobLength-shotNum(i)*postShootCut*Fs)*100;
    if partialSatis(i) < 0
        partialSatis(i) = 0;
    end
end

partialSatis

count = 0;
missionLength = imuJobEndIdx(end) - syncLen - imuJobStartIdx(1);
for i = imuJobStartIdx(1):imuJobEndIdx(end)-syncLen
    if (~isempty(find(totalVelocity(i:i+syncLen)>threshold)))
        if isempty(find((photoTimeS - imuTimeS(i)) .* (cutEnd - imuTimeS(i))<0))
            count = count + 1;
        end
    end
end
totalSatis = count / (missionLength-syncLen-sum(shotNum)*postShootCut*Fs) * 100;

totalSatis

partialVel = []; %%%TODO partial vel acquirement
partialVelBar = [];
totalVelMem = 0;
totalJobLength = 0;
for i = 1:length(imuJobStartIdx)
    jobLength = 0;
    velMem = 0;
    velBarMem = zeros(size(velBar,1)-1,1);
    for j = imuJobStartIdx(i):imuJobEndIdx(i)
        if isempty(find((photoTimeS - imuTimeS(j)) .* (cutEnd - imuTimeS(j))<0))
            jobLength = jobLength + 1;
            totalJobLength = totalJobLength + 1;
            velMem = velMem + totalVelocity(j);
            velBarMem = velBarMem + velBar(1:end-1,j);
            totalVelMem = totalVelMem + totalVelocity(j);
        end
    end
    partialVelBar(:,i) = horzcat(velBarMem/jobLength);
    partialVel(i) = velMem / jobLength;
end

totalVel = totalVelMem / totalJobLength
partialVel

% save('partialVel_sortie2','partialVel');

%% job satisfactory rate trend

partialSatis1 = load('partialSatis_sortie1.mat');
partialSatis1 = partialSatis1.partialSatis;
partialSatis2 = load('partialSatis_sortie2.mat');
partialSatis2 = partialSatis2.partialSatis;
partialVel1 = load('partialVel_sortie1.mat');
partialVel1 = partialVel1.partialVel;
partialVel2 = load('partialVel_sortie2.mat');
partialVel2 = partialVel2.partialVel;

shotNum1 = [0 20 0 18 0 21 0 27 0 19 0 20];
shotNum2 = [0 20 0 19 0 22 0 27 0 20 0 19];

photoCount1 = shotNum1;
faultPhotoCount1 = [0 4 0 1 0 1 0 16 0 3 0 3];
realBadPhotoCount1 = [0 1 0 1 0 0 0 7 0 0 0 1];

photoCount2 = shotNum2;
faultPhotoCount2 = [0 0 0 1 0 1 0 3 0 1 0 3];
realBadphotoCount2 = [0 0 0 0 0 0 0 1 0 0 0 0];

faultRatio1 = faultPhotoCount1 ./ photoCount1 * 100;
realBadRatio1 = realBadPhotoCount1 ./ photoCount1 * 100;
faultRatio2 = faultPhotoCount2 ./ photoCount2 * 100;
realBadRatio2 = realBadphotoCount2 ./ photoCount2 * 100;

figure(12)
clf
hold on
grid on
scatter(horzcat(partialSatis1,partialSatis2),horzcat(faultRatio1,faultRatio2),'LineWidth',5)
scatter(horzcat(partialSatis1,partialSatis2),horzcat(realBadRatio1,realBadRatio2),'LineWidth',5)
plot([0 60], [0 60], 'k')
scatter([9.3789 2.197],[22.4 7.08],'x','LineWidth',4)
% axis equal
xlim([0 60])
ylim([0 60])

disp('real violation rate')
sum(faultPhotoCount1)/sum(photoCount1)*100
sum(faultPhotoCount2)/sum(photoCount2)*100

xlabel('Expected fault ratio [%]')
ylabel('Actual fault ratio [%]')
title('Fault ratio expectation - actual comparison')

figure(1212)
clf
hold on
grid on
scatter(horzcat(partialVel1,partialVel2),horzcat(faultRatio1,faultRatio2),'LineWidth',5)
scatter(horzcat(partialVel1,partialVel2),horzcat(realBadRatio1,realBadRatio2),'LineWidth',5)
scatter([1.697 1.2517],[22.4 7.08],'x','LineWidth',4)
% plot([0 60], [0 60], 'k')
% axis equal
% xlim([0 60])
% ylim([0 60])

xlabel('avg pt vel [%]')
ylabel('Actual fault ratio [%]')
title('Fault avg ptvel - actual comparison')

%% photo failure comparison
preShootWindow = 0.1;
shootAvgTotVal = [];

for i = 1:length(photoTimeS)
    shootTimeIdx = findNearestIdx(imuTimeS,photoTimeS(i));
    preShootWindowIdx = findNearestIdx(imuTimeS,photoTimeS(i) - preShootWindow);
    shootAvgTotVal(i) = mean(totalVelocity(preShootWindowIdx:shootTimeIdx));
end

figure(14)
clf
hold on
grid on
photoScore = ones(1,length(photoTime));
photoScore(faultPhotoIdx) = 0.5;
photoScore(realBadPhotoIdx) = 0;
scatter(shootAvgTotVal, photoScore)

%% comparison plot
figure(888)
subplot(3,1,3)
clf
hold on
grid on
plot(imuTimeS,totalVelocity)

kernel = [2 4 6 8 10 12];
figure(13)
clf
bar([partialSatis1(kernel)' , partialSatis2(kernel)'])
hold on
plot(faultRatio1(kernel),'bx','MarkerSize',10,'LineWidth',4)
plot(faultRatio2(kernel),'rx','MarkerSize',10,'LineWidth',4)
grid on
legend('sortie 1 (ouster 20hz)','sortie 2 (ouster 10hz)')
xlabel('job #')
ylabel('expected fault ratio')
