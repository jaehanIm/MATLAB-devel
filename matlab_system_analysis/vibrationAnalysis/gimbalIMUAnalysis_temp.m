addpath('./..')
%% Load data


% 220210 ss1
imuData = readtable('/home/jaehan/log/220210_temp/220210_141454_soso1/aSensorImu_220210_141454.csv');
imuData(1:165,:) = [];
gdLog = readtable('/home/jaehan/log/220210_temp/220210_141454_soso1/gdLog_220210_141454.csv');
photoDir = dir('/home/jaehan/log/220210_temp/220210_141454_soso1/soft-soft sortie 1');
photoXml = readtable('/home/jaehan/log/220210_temp/220210_141454_soso1/df_xmp.csv');
photoTimeTrim = seconds(0.8 + 0.6780 - 0.1);
% imuData(1680000:end,:) = []; % 6-1
imuData(1:1680000,:) = []; % 6-2
% gdLog(21600:end,:) = []; % 6-1
gdLog(1:21600,:) = []; % 6-2

% 220210 st1
% imuData = readtable('/home/jaehan/log/220210_temp/220210_150409_stso1/aSensorImu_220210_150409.csv');
% imuData(1:165,:) = [];
% gdLog = readtable('/home/jaehan/log/220210_temp/220210_150409_stso1/gdLog_220210_150409.csv');
% photoDir = dir('/home/jaehan/log/220210_temp/220210_150409_stso1/stiff-soft sortie 1');
% photoXml = readtable('/home/jaehan/log/220210_temp/220210_150409_stso1/df_xmp.csv');
% photoTimeTrim = seconds(0.8 + 0.6780 - 0.1);
% % imuData(1731564:end,:) = []; % 6-1
% imuData(1:1731564,:) = []; % 6-2
% % gdLog(19179:end,:) = [];
% gdLog(1:19179,:) = [];

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

photoJobCategory = [];
for i = 1:length(jobStartIdx)
    photoJobCategory{i} = find(photoTimeS > jobStartIdxS(i) & photoTimeS < jobEndIdxS(i));
end

disp("Data Loading Complete!")

%% pre processing
acc_0 = detrend(acc_0,1);
acc_1 = detrend(acc_1,1);
acc_2 = detrend(acc_2,2);
gyro_0 = detrend(gyro_0,1);
gyro_1 = detrend(gyro_1,1);
gyro_2 = detrend(gyro_2,1);

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
% for i = 1:length(jobStartIdx)
%     text(jobStartIdxS(i)+2,3,jobName(i),'Fontsize',10,'FontWeight','bold','Color','m')
% end
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

% for i = 1:length(jobStartIdx)
%     text(jobStartIdxS(i)+2,30,jobName(i),'Fontsize',10,'FontWeight','bold','Color','m')
% end

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

totalVel = totalVelMem / totalJobLength;
partialVel

%% individual level
preShootWindow = 0.1;
shootAvgTotVal = [];

for i = 1:length(photoTimeS)
    shootTimeIdx = findNearestIdx(imuTimeS,photoTimeS(i));
    preShootWindowIdx = findNearestIdx(imuTimeS,photoTimeS(i) - preShootWindow);
    shootAvgTotVal(i) = mean(totalVelocity(preShootWindowIdx:shootTimeIdx));
end

% save('partialVel_sortie_ss','partialVel');
% save('partialSatis_sortie_ss','partialSatis');
% save('totalSatis_sortie_ss','totalSatis');
% save('totalVel_sortie_ss','totalVel');