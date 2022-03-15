addpath('./..')
%% yongdaeri
sample{1} = [2 4 12 16 24 28 30 32 40 67 92 105 109 110 111 112 116 121 122 131 138 140 169 174 178 182 192 194 204 206 210 230 232 234 244 248 261 263 265 268 270 272 274 281 283 286 288 305 307];
sample{2} = [2 6 12 16 32 121 122 169 174 175 178 180 182 184 200 203 230 232 258 259 260 261 265 266 268 269 270 272 274 281 283 286 290 295 310 311 314 316 318 320 324];
sample{3} = [2, 12, 16, 32, 121, 122, 169, 174, 178, 180, 182, 230, 281, 288, 295, 315, 328];

single = []; double = []; triple = [];

for i = 1:341
    count = 0;
    for j = 1:3
        if ismember(i,sample{j})
            count = count + 1;
        end
    end
    if count == 1
        single = vertcat(single,i);
    elseif count == 2
        double = vertcat(double,i);
    elseif count == 3
        triple = vertcat(triple,i);
    end
end

%% Load data

% 220210 ss1
gdLog1 = readtable('/home/jaehan/log/220210_temp/220210_141454_soso1/gdLog_220210_141454.csv');
photoXml = readtable('/home/jaehan/Desktop/yongdaeri/yongdaeri_xmp.csv');
photoTimeTrim = seconds(0.8 + 0.6780 - 0.1);

% 220210 st1
gdLog2 = readtable('/home/jaehan/log/220210_temp/220210_150409_stso1/gdLog_220210_150409.csv');
gdLog = [gdLog1;gdLog2];

%% Preprocessing

photoTime = photoXml.xmpTime;
photoTime = datetime(photoTime,'TimeZone','Asia/Tokyo');
photoTime = photoTime + photoTimeTrim;
photoName = photoXml.imgNumber;

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

% time synchronization
gdTime = datetime(gdLog.rosTime(1:end),'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
photoTimeDelay = photoTime(1) - gdTime(1);
photoTimeS = seconds(photoTime - photoTime(1) + photoTimeDelay);
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);
photoTimeS(photoTimeS<=0) = 0;

jobStartIdxS = gdTimeS(jobStartIdx);
jobEndIdxS = gdTimeS(jobEndIdx);


%% Categorization
photoJobCategory = [];
for i = 1:length(jobStartIdx)
    photoJobCategory{i} = find(photoTimeS > jobStartIdxS(i) & photoTimeS < jobEndIdxS(i));
end

singleShotTime = [];
doubleShotTime = [];
tripleShotTime = [];
for i = 1:length(single)
    singleShotTime(i) = photoTimeS(find(photoName == single(i)));
end
for i = 1:length(double)
    doubleShotTime(i) = photoTimeS(find(photoName == double(i)));
end
for i = 1:length(triple)
    tripleShotTime(i) = photoTimeS(find(photoName == triple(i)));
end

singlePhotoJobCategory = [];
doublePhotoJobCategory = [];
triplePhotoJobCategory = [];
for i = 1:length(jobStartIdxS)
    singlePhotoJobCategory{i} = find(singleShotTime > jobStartIdxS(i) & singleShotTime < jobEndIdxS(i));
    doublePhotoJobCategory{i} = find(doubleShotTime > jobStartIdxS(i) & doubleShotTime < jobEndIdxS(i));
    triplePhotoJobCategory{i} = find(tripleShotTime > jobStartIdxS(i) & tripleShotTime < jobEndIdxS(i));
end

singleCount = zeros(50,1); doubleCount = zeros(50,1); tripleCount = zeros(50,1);
for i = 1:length(jobStartIdx)
    singleCount(i) = size(singlePhotoJobCategory{i},2);
    doubleCount(i) = size(doublePhotoJobCategory{i},2);
    tripleCount(i) = size(triplePhotoJobCategory{i},2);
end

countScore = mean([singleCount, doubleCount, tripleCount],2);

disp("Data Loading Complete!")

partialVel_ss = load('partialVel_sortie_ss.mat'); partialVel_ss = partialVel_ss.partialVel;
partialVel_st = load('partialVel_sortie_st.mat'); partialVel_st = partialVel_st.partialVel;
partialSatis_ss = load('partialSatis_sortie_ss.mat'); partialSatis_ss = partialSatis_ss.partialSatis;
partialSatis_st = load('partialSatis_sortie_st.mat'); partialSatis_st = partialSatis_st.partialSatis;

partialVel = [partialVel_ss,partialVel_st];
partialSatis = [partialSatis_ss,partialSatis_st];

partialVel = partialVel';
partialSatis = partialSatis';

%% mission level
totalVel = [1.0008,1.1854,1.2438,1.3051];

missionScore= [mean(nonzeros(countScore(1:12))),mean(nonzeros(countScore(13:25))),mean(nonzeros((countScore(26:37)))),mean(nonzeros(countScore(38:50)))];

figure(3)
clf
hold on
grid on
plot(totalVel,missionScore,'o','LineWidth',2)
xlabel('mission average pt.vel.')
ylabel('mission average fault score')
title('HI fault score to pt.vel. (mission average)')

corrcoef(totalVel,missionScore)

%% shoot level
ind1 = load('shootAvgTotVal_ss.mat');
ind2 = load('shootAvgTotVal_st.mat');
ind1 = ind1.shootAvgTotVal;
ind2 = ind2.shootAvgTotVal;
ind = [ind1,ind2]';

% individual score
indCountScore = zeros(341,1);
for i = 1:341
    if ismember(i,single)
        indCountScore(i) = 1;
    elseif ismember(i,double)
        indCountScore(i) = 2;
    elseif ismember(i,triple)
        indCountScore(i) = 3;
    end
end

figure(4)
clf
hold on
grid on
% plot(ind(indCountScore ~= 0),nonzeros(indCountScore),'o','LineWidth',2)
plot(ind,indCountScore,'o','LineWidth',2)
xlabel('pt.vel.')
ylabel('individual fault score')
title('HI fault score to pt.vel. (individual shoot)')

corrcoef(ind,indCountScore)
%% job level

partialVel(countScore == 0) = [];
partialSatis(countScore == 0) = [];
countScore(countScore == 0) = [];

figure(1)
clf
hold on
grid on
plot(partialVel,countScore,'o','LineWidth',2)
xlabel('job average pt.vel.')
ylabel('job average fault score')
title('HI fault score to pt.vel. (job average)')

figure(2)
clf
hold on
grid on
plot(partialSatis,countScore,'o')

corrcoef(partialVel,countScore)
corrcoef(partialSatis,countScore)


