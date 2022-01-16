addpath('./..')

%% Load data

% 211223 oksang
% imuData = readtable('/home/jaehan/log/211223_[1]_oksang/aSensorLog_211223_120453.csv');
% gdLog = readtable('/home/jaehan/log/211223_[1]_oksang/gdLog_211223_120453.csv');
% photoXml = readtable('/home/jaehan/log/211223_[1]_oksang/df_xmp.csv');
% photoBlurEst = readtable('/home/jaehan/log/211223_[1]_oksang/blurTestResult.csv');
% photoTimeTrim = seconds(1.179);

% 220105 oksang
% imuData = readtable('/home/jaehan/log/220105_[1]_oksang/220105_133822/aSensorLog_220105_133822.csv');
% gdLog = readtable('/home/jaehan/log/220105_[1]_oksang/220105_133822/gdLog_220105_133822.csv');
% photoXml = readtable('/home/jaehan/log/220105_[1]_oksang/220105_133822/df_xmp.csv');
% opts = detectImportOptions('/home/jaehan/log/220105_[1]_oksang/220105_133822/s1_output.csv');
% opts.DataLines = 2;
% photoBlurEst = readtable('/home/jaehan/log/220105_[1]_oksang/220105_133822/s1_output.csv',opts);
% photoTimeTrim = seconds(1.179-0.626-0.1);

% 220107 oksang - 20hz
imuData = readtable('/home/jaehan/log/220107_[1]_oksang/220107_152745_20hz/aSensorLog_220107_152745.csv');
gdLog = readtable('/home/jaehan/log/220107_[1]_oksang/220107_152745_20hz/gdLog_220107_152745.csv');
photoXml = readtable('/home/jaehan/log/220107_[1]_oksang/220107_152745_20hz/df_xmp.csv');
opts = detectImportOptions('/home/jaehan/log/220107_[1]_oksang/220107_152745_20hz/output_2.csv');
opts.DataLines = 2;
photoBlurEst = readtable('/home/jaehan/log/220107_[1]_oksang/220107_152745_20hz/output_2.csv',opts);
photoTimeTrim = seconds(1.179-0.626+0.366-0.1);
manPhotoFault1 = [014  015  016  017  018 022  023  024  027  028 029  031  032  033  034 036  038  040  065  069 076  080  083  084  086 087  089  090  091  092 093  094  096  097  099 100  132  135  161  165 169  177  180  189  191 193  196  199]- 4;
manPhotoFault2 = [13 15 16 17 18 22 23 24 27 28 29 30 32 33 34 36 38 40 48 65 66 69 72 76 80 83 84 86 87 89 90 91 92 93 94 96 97 99 100 103 112 116 132 133 135 136 142 161 164 165 169 177 180 189 191 193 195 196 199]- 4;
manPhotoFault3 = [14, 15, 16, 17, 18, 22, 23, 24, 27, 28, 29, 31, 32, 33, 34, 38, 66, 69, 76, 80, 83, 84, 86, 87, 89, 90, 91, 92, 93, 94, 96, 97, 99, 100, 132, 135, 161, 165, 169, 177, 180, 189, 191, 193, 196, 199]- 4;
manPhotoFaultTotal = union(union(manPhotoFault1,manPhotoFault2),manPhotoFault3);
singleFault = []; doubleFault = []; tripleFault = [];
for i = 1:length(manPhotoFaultTotal)
    candidate = manPhotoFaultTotal(i);
    detectNum = sum([ismember(candidate,manPhotoFault1),ismember(candidate,manPhotoFault2),ismember(candidate,manPhotoFault3)]);
    if detectNum == 1
        singleFault = horzcat(singleFault,candidate);
    elseif detectNum == 2
        doubleFault = horzcat(doubleFault,candidate);
    else
        tripleFault = horzcat(tripleFault,candidate);
    end
end

% 220107 oksang - 10hz
% imuData = readtable('/home/jaehan/log/220107_[1]_oksang/220107_150302_10hz/aSensorLog_220107_150302.csv');
% gdLog = readtable('/home/jaehan/log/220107_[1]_oksang/220107_150302_10hz/gdLog_220107_150302.csv');
% photoXml = readtable('/home/jaehan/log/220107_[1]_oksang/220107_150302_10hz/df_xmp.csv');
% opts = detectImportOptions('/home/jaehan/log/220107_[1]_oksang/220107_150302_10hz/output_1.csv');
% opts.DataLines = 2;
% photoBlurEst = readtable('/home/jaehan/log/220107_[1]_oksang/220107_150302_10hz/output_1.csv');
% photoTimeTrim = seconds(1.179-0.626+0.39);
% manPhotoFault1 = [003  009  014  015  016 018  020  021  026  027 029  032  037  042  043 050  057  071  072  073 075  079  082  083  085 090  093  095  096  097 101  149  155  157  159 164  170  176  184]- 2;
% manPhotoFault2 = [3 9 10 11 14 15 16 18 20 21 23 27 29 32 33 35 37 42 43 50 57 70 71 72 73 75 77 80 83 85 90 93 95 96 101 106 109 155 157 164 170 173 180 186 197]- 2;
% manPhotoFault3 = [3 9 14 15 16 20 21 26 27 29 32 42 43 59 63 71 72 73 75 77 79 83 85 90 93 95 97 101 149 155 159 164 170 173 184 196]- 2;
% manPhotoFaultTotal = union(union(manPhotoFault1,manPhotoFault2),manPhotoFault3);
% singleFault = []; doubleFault = []; tripleFault = [];
% for i = 1:length(manPhotoFaultTotal)
%     candidate = manPhotoFaultTotal(i);
%     detectNum = sum([ismember(candidate,manPhotoFault1),ismember(candidate,manPhotoFault2),ismember(candidate,manPhotoFault3)]);
%     if detectNum == 1
%         singleFault = horzcat(singleFault,candidate);
%     elseif detectNum == 2
%         doubleFault = horzcat(doubleFault,candidate);
%     else
%         tripleFault = horzcat(tripleFault,candidate);
%     end
% end

% photoMaxVal?
for i = 1:size(photoBlurEst,1)
    if ~isempty(photoBlurEst.totalRecord{i})
        tempRecord = str2num(photoBlurEst.totalRecord{i});
        photoBlurEst.maxVal(i) = max(tempRecord);
        photoBlurEst.meanVal(i) = mean(tempRecord);
    else
        photoBlurEst.maxVal(i) = nan;
        photoBlurEst.maxVal(i) = nan;
        photoBlurEst.meanVal(i) = nan;
    end
end

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
    elseif (gdLog.fcMcMode(i) == 1 && gdLog.fcMcMode(i-1) == 2) ||(gdLog.fcMcMode(i) == 0 && gdLog.fcMcMode(i-1) == 2)
        jobEndIdx = horzcat(jobEndIdx,i);
    elseif gdLog.fcMcMode(i) == 255 && gdLog.fcMcMode(i-1) == 0
        landIdx = i;
    end
end

% TEMPORARY jobEndidx
% jobEndIdx = [jobEndIdx,size(gdLog,1)];


% jobName = ["Sweep","FB", "LR", "UD", "YawRL", "Hover"];
% jobName = ["FB", "LR", "UD", "YawRL", "Hover"];
jobName = ["Sweep_FM","Sweep_AM","FB", "LR", "UD", "YawRL", "Hover"];
parseStart = 1;
parseEnd = size(imuData,1);

% time synchronization
imuTime = datetime(imuData.rosTime(parseStart:parseEnd),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
gdTime = datetime(gdLog.rosTime(1:end),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
imuTimeDelay = seconds(imuTime(1) - gdTime(1));
photoTimeDelay = photoTime(1) - gdTime(1);
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);
imuTimeS = seconds(imuTime - imuTime(1)) + imuTimeDelay;
photoTimeS = seconds(photoTime - photoTime(1) + photoTimeDelay);
photoTimeS(photoTimeS<=0) = 0;

% TEMPORARY time synchronization
% imuTime = datetime(imuData.rosTime,'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
% gdTime = datetime(gdLog.rosTime,'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
% timeError = imuTime(1) - gdTime(1);
% imuTime = imuTime - timeError;
% timeDelay = 31.9;
% photoTimeDelay = 0.9;
% gdTimeS = (gdLog.rosTime - gdLog.rosTime(1));
% imuTimeS = seconds(imuTime - seconds(timeDelay) - imuTime(1));
% photoTimeS = seconds(photoTime - gdTime(1) + seconds(photoTimeDelay));
% photoTimeS(photoTimeS<=0) = 0;



jobStartIdxS = gdTimeS(jobStartIdx);
jobEndIdxS = gdTimeS(jobEndIdx);

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
    interest = find(imuTimeS>=gdTimeS(i) & imuTimeS<=gdTimeS(i+1));
    
    gimbalRF = gdLog.gimbalRpy_deg_0(i);
    gimbalRB = gdLog.gimbalRpy_deg_0(i+1);
    gimbalPF = gdLog.gimbalRpy_deg_1(i);
    gimbalPB = gdLog.gimbalRpy_deg_1(i+1);
    gimbalY = 0;
    dcmI2GF = angle2dcm(gimbalY, gimbalRF, gimbalPF,'ZXY');
    dcmI2GB = angle2dcm(gimbalY, gimbalRB, gimbalPB,'ZXY');
    
    timeF = gdTimeS(i);
    timeB = gdTimeS(i+1);
    gduF = gdLog.velUVW_mps_0(i);
    gduB = gdLog.velUVW_mps_0(i+1);
    gdvF = gdLog.velUVW_mps_1(i);
    gdvB = gdLog.velUVW_mps_1(i+1);
    gdwF = gdLog.velUVW_mps_2(i);
    gdwB = gdLog.velUVW_mps_2(i+1);
    
    gdVelF = dcmI2GF * [gduF;gdvF;gdwF];
    gdVelB = dcmI2GB * [gduB;gdvB;gdwB];
    
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
for i = 1:length(jobStartIdxS)
plot([jobStartIdxS(i) jobStartIdxS(i)],[-std(v)*5 std(v)*5],'r--')
plot([jobEndIdxS(i) jobEndIdxS(i)],[-std(v)*5 std(v)*5],'g--')
end
legend('integrated v (HP filtered)','gdLog v')

disp('Partial FDI complete!')
%% data selection
y = acc_0;

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
% figure(7)
% clf
% hold on
% grid on
% Y = fft(y);
% 
% freqL = [0.01 5  30  100 500 800 1200];
% freqH = [5 30 100 500 800 1200 2000];
% df = L/Fs;
% filteredFreq_R = [];
% 
% for i = 1:length(freqL)
%     P_filter = Y;
%     kernel = zeros(L,1);
%     kernel(freqL(i)*df:freqH(i)*df) = 1;
%     kernel(L-freqH(i)*df:L-freqL(i)*df) = 1;
%     P_filter = P_filter .* kernel;
%     fftRegen = ifft(P_filter);
%     
%     filteredFreq_R(i,:) = real(fftRegen);
%     
%     subplot(length(freqL),1,i)
%     plot(imuTimeS,real(fftRegen))
% %     xlim([10 11])
%     labelText = num2str(freqL(i)) + "hz to " + num2str(freqH(i)) + "hz";
%     title(labelText);
% end
% 
% sgtitle('FFT spectral decompostion','fontsize',14) 

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
for i = 1:length(jobStartIdxS)
    plot([0 200],[jobStartIdxS(i) jobStartIdxS(i)],'m--','LineWidth',1.5)
    plot([0 200],[jobEndIdxS(i) jobEndIdxS(i)],'g--','LineWidth',1.5)
    text(2,jobStartIdxS(i)+2,jobName(i),'Color',[1 1 1])
end


subplot(1,4,2)
imagesc([40 500],[imuTimeS(1) imuTimeS(end)],spectrogram(:,40*timeStep:500*timeStep));
hold on
for i = 1:length(jobStartIdxS)
    plot([40 500],[jobStartIdxS(i) jobStartIdxS(i)],'m--','LineWidth',1.5)
    plot([40 500],[jobEndIdxS(i) jobEndIdxS(i)],'g--','LineWidth',1.5)
end
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
for i = 1:length(jobStartIdxS)
    plot([jobStartIdxS(i) jobStartIdxS(i)],[-100 100],'m--','LineWidth',1.5)
    plot([jobEndIdxS(i) jobEndIdxS(i)],[-100 100],'g--','LineWidth',1.5)
    text(jobStartIdxS(i)+2,2,jobName(i),'Color',[1 1 1])
end

disp("Spectrogram plotting complete!")


%% Vibration component analysis - pqr
p = gyro_0;
q = gyro_1;
r = gyro_2;

P = fft(p);
Q = fft(q);
R = fft(r);

% freqL = [0 5  30  100 500 800 1200];
% freqH = [5 30 100 500 800 1200 2000];
freqL = [0 5  30  100 500 1200];
freqH = [5 30 100 500 800 2000];
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

% freqL = [0 5  30  100 500 800 1200];
% freqH = [5 30 100 500 800 1200 2000];
freqL = [0 5  30  100 500 1200];
freqH = [5 30 100 500 800 2000];
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


vel = [];
x_comp = [];
y_comp = [];
pixelPosAngle = pi/4;
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
for i = 1:size(jobStartIdx,2)
    plot([jobStartIdxS(i) jobStartIdxS(i)], [0 max(totalVelocity)],'r--')
    plot([jobEndIdxS(i) jobEndIdxS(i)], [0 max(totalVelocity)],'g--')
end

for i = 1:length(photoTimeS)
    plot(photoTimeS(i),0,'g*','LineWidth',4,'MarkerSize',10)
end

for i = 1:length(jobStartIdx)
    text(jobStartIdxS(i)+2,30,jobName(i),'Fontsize',10,'FontWeight','bold','Color','m')
end

plot(photoTimeS,photoBlurEst.minVal,'rx--','LineWidth',3,'MarkerSize',8)

legend('pointing velocity','threshold');
ylim([0 40])


%% Total contribution

% start = imuJobStartIdx(7);
% finish = imuJobEndIdx(7);

start = imuJobEndIdx(6);
finish = L;

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
(count- sum(shotNum) * 560)/(L-20 - sum(shotNum) * 560) * 100


%% Partial satisfaction rate
postShootCut = 0.34;
cutEnd = photoTimeS + postShootCut;
syncLen = 3;

partialSatis = [];
for i = 1:length(imuJobStartIdx)
    jobLength = imuJobEndIdx(i) - imuJobStartIdx(i);
    count = 0;
    for j = imuJobStartIdx(i):imuJobEndIdx(i)-syncLen
        if (~isempty(find(totalVelocity(j:j+syncLen)>threshold)))
            if isempty(find((photoTimeS - imuTimeS(j)) .* (cutEnd - imuTimeS(j))<0))
                count = count + 1;
            end
        end
    end
%     partialSatis(i) = (count - shotNum(i) * 560)/(jobLength-20-shotNum(i) * 560)*100;
    partialSatis(i) = count/(jobLength-syncLen-shotNum(i)*postShootCut*Fs)*100;
    if partialSatis(i) < 0
        partialSatis(i) = 0;
    end
end

partialSatis

partialVel = []; %%%TODO partial vel acquirement
partialVelBar = [];
for i = 1:length(imuJobStartIdx)
    jobLength = 0;
    velMem = 0;
    velBarMem = zeros(size(velBar,1)-1,1);
    for j = imuJobStartIdx(i):imuJobEndIdx(i)
        if isempty(find((photoTimeS - imuTimeS(j)) .* (cutEnd - imuTimeS(j))<0))
            jobLength = jobLength + 1;
            velMem = velMem + totalVelocity(j);
            velBarMem = velBarMem + velBar(1:end-1,j);
        end
    end
    partialVelBar(:,i) = horzcat(velBarMem/jobLength);
    partialVel(i) = velMem / jobLength;
end

partialVel

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

figure(19)
bar((partialVelBar)','stacked')
title('contribution portion')
xlabel('job #')
ylabel('m/s')
figure(1919)
bar((partialVelBar./sum(partialVelBar)*100)','stacked')
title('normalized contribution portion')
xlabel('job #')
ylabel('m/s')

% temp = partialVelBar(:,6);
% save('20hz_partialVel2.mat','temp');

%% Photo shoot time considered trend analysis

preShootWindow = 0.1;
shootAvgTotVal = [];

for i = 1:length(photoTimeS)
    shootTimeIdx = findNearestIdx(imuTimeS,photoTimeS(i));
    preShootWindowIdx = findNearestIdx(imuTimeS,photoTimeS(i) - preShootWindow);
    shootAvgTotVal(i) = mean(totalVelocity(preShootWindowIdx:shootTimeIdx));
end

photoScoreMin = photoBlurEst.minVal;
photoScoreMax = photoBlurEst.maxVal;
photoScoreMean = photoBlurEst.meanVal;
% for i = 1:length(photoScoreMin)
%     if isnan(photoScoreMin(i))
%         photoScoreMin(i) = 0;
%         photoScoreMax(i) = 0;
%         photoScoreMean(i) = 0;
%     end
% end 

% omit index
omitIdx = find(isnan(photoScoreMin));


figure(12)
clf
hold on
grid on
plot(shootAvgTotVal,photoScoreMax,'*');
xlim([shootAvgTotVal(1)+0.1 inf])

partialFault_sync = [];
partialFault_dirac = [];
for j = 1:length(photoTimeS)
    count = 0;
    count2 = 0;
    shootTimeIdx = findNearestIdx(imuTimeS,photoTimeS(j));
    preShootWindowIdx = findNearestIdx(imuTimeS,photoTimeS(j) - preShootWindow);
    windowLength = shootTimeIdx - preShootWindowIdx;
    for i= preShootWindowIdx:shootTimeIdx - syncLen
        if(~isempty(find(totalVelocity(i:i+syncLen) > threshold)))
            count = count + 1;
        end
    end
    for i = preShootWindowIdx:shootTimeIdx
        if totalVelocity(i) > threshold
            count2 = count2 + 1;
        end
    end
    partialFault_sync(j) = count / (windowLength - syncLen) * 100;
    partialFault_dirac(j) = count2 / windowLength * 100;
end

figure(13)
clf
hold on
grid on
plot(normalize(shootAvgTotVal,'range')*100,photoScoreMax,'*','MarkerSize',10,'LineWidth',3)
% plot(shootAvgTotVal,photoScoreMax,'*','MarkerSize',10,'LineWidth',3)
plot(partialFault_sync,photoScoreMax,'*','MarkerSize',10,'LineWidth',3)
plot(partialFault_dirac,photoScoreMax,'*','MarkerSize',10,'LineWidth',3)
legend('average tot val','sync fault','dirac fault')


figure(1212)
clf
hold on
grid on
plot(shootAvgTotVal)
plot(partialFault_sync)
plot(partialFault_dirac)
plot(photoScoreMin)

% job level trend analysis

jobPhotoAvgScoreMin = [];
jobPhotoAvgScoreMax = [];
jobPhotoAvgScoreMean = [];
jobPhotoStdMin = [];
jobPhotoStdMax = [];
jobPhotoStdMean = [];
for i = 1:length(photoJobCategory)
    jobPhotoAvgScoreMin(i) = mean(photoScoreMin(photoJobCategory{i}(1):photoJobCategory{i}(end)),'omitnan');
    jobPhotoAvgScoreMax(i) = mean(photoScoreMax(photoJobCategory{i}(1):photoJobCategory{i}(end)),'omitnan');
    jobPhotoAvgScoreMean(i) = mean(photoScoreMean(photoJobCategory{i}(1):photoJobCategory{i}(end)),'omitnan');
    
    jobPhotoStdMin(i) = std(photoScoreMin(photoJobCategory{i}(1):photoJobCategory{i}(end)),'omitnan');
    jobPhotoStdMax(i) = std(photoScoreMax(photoJobCategory{i}(1):photoJobCategory{i}(end)),'omitnan');
    jobPhotoStdMean(i) = std(photoScoreMean(photoJobCategory{i}(1):photoJobCategory{i}(end)),'omitnan');
end


figure(14)
clf
hold on
grid on
plot(partialSatis,jobPhotoAvgScoreMin,'ro','MarkerSize',10,'LineWidth',3);
plot(partialSatis,jobPhotoAvgScoreMean,'go','MarkerSize',10,'LineWidth',3);
plot(partialSatis,jobPhotoAvgScoreMax,'bo','MarkerSize',10,'LineWidth',3);

errorbar(partialSatis,jobPhotoAvgScoreMin,jobPhotoStdMin,'ro','MarkerSize',10,'LineWidth',1);
errorbar(partialSatis,jobPhotoAvgScoreMean,jobPhotoStdMean,'go','MarkerSize',10,'LineWidth',1);
errorbar(partialSatis,jobPhotoAvgScoreMax,jobPhotoStdMax,'bo','MarkerSize',10,'LineWidth',1);
plot(partialSatis(1),jobPhotoAvgScoreMin(1),'rx','MarkerSize',10,'LineWidth',3);
legend('Minimum','Average','Maximum')

title('Expected fault ratio to photo score','fontsize',14)
xlabel('EFR')
ylabel('photo score')

text(partialSatis+0.2,jobPhotoAvgScoreMean,jobName(1:length(partialSatis)),'Fontsize',10,'FontWeight','bold','Color','m')

%%% TODO
figure(1414)
clf
hold on
grid on
plot(partialVel,jobPhotoAvgScoreMin,'ro','MarkerSize',10,'LineWidth',3);
plot(partialVel,jobPhotoAvgScoreMean,'go','MarkerSize',10,'LineWidth',3);
plot(partialVel,jobPhotoAvgScoreMax,'bo','MarkerSize',10,'LineWidth',3);

errorbar(partialVel,jobPhotoAvgScoreMin,jobPhotoStdMin,'ro','MarkerSize',10,'LineWidth',1);
errorbar(partialVel,jobPhotoAvgScoreMean,jobPhotoStdMean,'go','MarkerSize',10,'LineWidth',1);
errorbar(partialVel,jobPhotoAvgScoreMax,jobPhotoStdMax,'bo','MarkerSize',10,'LineWidth',1);
plot(partialVel(1),jobPhotoAvgScoreMin(1),'rx','MarkerSize',10,'LineWidth',3);
legend('Minimum','Average','Maximum')

title('partial Velocity to photo score','fontsize',14)
xlabel('partial totVel')
ylabel('photo score')

text(partialVel+0.05,jobPhotoAvgScoreMax,jobName(1:length(partialVel)),'Fontsize',10,'FontWeight','bold','Color','m')



figure(15)
clf
hold on
grid on
plot(partialFault_sync,photoScoreMin,'*','MarkerSize',10,'LineWidth',3)
plot(partialFault_sync,photoScoreMean,'*','MarkerSize',10,'LineWidth',3)
plot(partialFault_sync,photoScoreMax,'*','MarkerSize',10,'LineWidth',3)
legend('min','mean','max')
title('partial fault sync base')

figure(1515)
clf
hold on
grid on
scoreTemp = shootAvgTotVal;
minTemp = photoScoreMin;
meanTemp = photoScoreMean;
maxTemp = photoScoreMax;
scoreTemp(photoJobCategory{5}) = [];
minTemp(photoJobCategory{5}) = [];
meanTemp(photoJobCategory{5}) = [];
maxTemp(photoJobCategory{5}) = [];
% plot(shootAvgTotVal,photoScoreMin,'*','MarkerSize',10,'LineWidth',3)
% plot(shootAvgTotVal,photoScoreMean,'*','MarkerSize',10,'LineWidth',3)
% plot(shootAvgTotVal,photoScoreMax,'*','MarkerSize',10,'LineWidth',3)
plot(scoreTemp, minTemp,'*','MarkerSize',10,'LineWidth',3)
plot(scoreTemp, meanTemp,'*','MarkerSize',10,'LineWidth',3)
plot(scoreTemp, maxTemp,'*','MarkerSize',10,'LineWidth',3)
legend('min','mean','max')
title('avg.pt.vel. to FFTS')
xlabel('avg.pt.vel. [m/s]')
ylabel('FFTS')

figure(151515)
clf
hold on
grid on
plot(partialFault_dirac,photoScoreMin,'*','MarkerSize',10,'LineWidth',3)
plot(partialFault_dirac,photoScoreMean,'*','MarkerSize',10,'LineWidth',3)
plot(partialFault_dirac,photoScoreMax,'*','MarkerSize',10,'LineWidth',3)
legend('min','mean','max')
title('partial fault dirac base')

figure(15151515)
clf
hold on
grid on
for i = 1:length(jobName)
    plot(shootAvgTotVal(photoJobCategory{i}),photoScoreMax(photoJobCategory{i}),'*','MarkerSize',10,'LineWidth',3)
end
legend(jobName(1),jobName(2),jobName(3),jobName(4),jobName(5),jobName(6),jobName(7))
title('avg.pt.vel. to FFTS')
xlabel('avg.pt.vel. [m/s]')
ylabel('FFTS')


% delete specific job
% for i = 1:length(jobName)
%     photoJobCategory{i}(find(ismember(photoJobCategory{i},omitIdx))) = [];
% end

% TODO
disp('partial corrcoef')
for i = 1:length(jobName)
    corrcoef(shootAvgTotVal(photoJobCategory{i}),photoScoreMax(photoJobCategory{i}))
end

% disp('for back')
% interest = [photoJobCategory{2};photoJobCategory{3};photoJobCategory{4};photoJobCategory{6};photoJobCategory{7}];
% maxCorr = corrcoef([partialFault_sync(interest)',photoScoreMax(interest)])


figure(16)
clf
% subplot(1,3,1)
hold on
grid on
manPhotoScore = zeros(1,length(shootAvgTotVal));
manPhotoScore(singleFault) = 1;
manPhotoScore(doubleFault) = 2;
manPhotoScore(tripleFault) = 3;
% plot(shootAvgTotVal(1:196),manPhotoScore(1:196),'o','MarkerSize',6,'LineWidth',3)
% xlabel('average velocity')

% subplot(1,3,2)
% hold on
% grid on
% plot(partialFault_sync(1:196),manPhotoScore(1:196),'o','MarkerSize',10,'LineWidth',3)
% xlabel('partial Fault ratio (sync)')

% subplot(1,3,3)
hold on
grid on
plot(photoScoreMax(1:196),manPhotoScore(1:196),'o','MarkerSize',10,'LineWidth',3)
xlabel('fft score')


partialFault_sync(omitIdx) = []; photoScoreMin(omitIdx) = []; photoScoreMax(omitIdx) = []; photoScoreMean(omitIdx) = []; shootAvgTotVal(omitIdx) = [];
% 
% 
% disp('sync fault')
% maxCorr = corrcoef([partialFault_sync',photoScoreMax])
% MeanCorr = corrcoef([partialFault_sync',photoScoreMean])
% minCorr = corrcoef([partialFault_sync',photoScoreMin])
% 
% disp('average tot val')
% maxCorr = corrcoef([shootAvgTotVal',photoScoreMax])
% MeanCorr = corrcoef([shootAvgTotVal',photoScoreMean])
% minCorr = corrcoef([shootAvgTotVal',photoScoreMin])
% 
% disp('man fault corr')
% a = corrcoef([shootAvgTotVal(1:196)',manPhotoScore(1:196)'])
% a = corrcoef([partialFault_sync(1:196)',manPhotoScore(1:196)'])
% a = corrcoef([photoScoreMax(1:196),manPhotoScore(1:196)'])
