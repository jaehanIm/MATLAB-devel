% function MC_analysis()
addpath('./data/')

%% Load data
dat = load('set4.mat');
dat = dat.MCData;
% dat(end,:,:) = [];
%% Analysis

TestClusteringTime = zeros(size(dat));
TestCompleteTime = zeros(size(dat));
TestSolveTime = zeros(size(dat));
TestScore = zeros(size(dat));
TestScoreL = zeros(size(dat));
CompCompleteTime = zeros(size(dat));
CompSolveTime = zeros(size(dat));
CompScore = zeros(size(dat));
CompScoreL = zeros(size(dat));
degreeConn = zeros(size(dat));
N = zeros(size(dat));
equiPerformanceTime = zeros(size(dat));
equiPerformanceTimeL = zeros(size(dat));
equiTimePerformance = zeros(size(dat));

iterN = size(dat,3);

% fovFactorSet = [3:-0.5:1,0.75];
% conThresSet = 5:10:55;
% sizeList = 3:-0.5:0.5;
% connList = 5:5:50;
% [X,Y] = meshgrid(conThresSet,fovFactorSet);

% parsing
for sizeFactor = 1:size(dat,1)
    for connFactor = 1:size(dat,2)
        for i = 1:iterN
            TestClusteringTime(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.TestClusteringTime;
            TestCompleteTime(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.TestInterCompleteTime;
            TestSolveTime(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.TestSolveTime;
            TestScore(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.TestScoreHist;
            TestScoreL(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.TestScoreHistL;
            CompCompleteTime(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.CompCompleteTime;
            CompSolveTime(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.CompSolveTime;
            CompScore(sizeFactor,connFactor,i) = max(dat{sizeFactor,connFactor,i}.CompScorePer);
            CompScoreL(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.CompScoreL;
            degreeConn(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.degreeConnectivity;
            N(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.N;
            if ~isempty(dat{sizeFactor,connFactor,i}.equiPerformanceTimeL)
%                 equiPerformanceTime(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.equiPerformanceTime;
                equiPerformanceTimeL(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.equiPerformanceTimeL;
            else
%                 equiPerformanceTime(sizeFactor,connFactor,i) = -1;
                equiPerformanceTimeL(sizeFactor,connFactor,i) = 0;
            end
            equiTimePerformance(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.equiTimePerformance;
        end
    end
end

TestTotalTime = TestClusteringTime + TestCompleteTime + TestSolveTime;
CompTotalTime = CompCompleteTime + CompSolveTime;


% mean - stdev
TestClusteringTime = meanStdProcess(TestClusteringTime);
TestCompleteTime = meanStdProcess(TestCompleteTime);
TestSolveTime = meanStdProcess(TestSolveTime);
TestScore = meanStdProcess(TestScore);
TestScoreL = meanStdProcess(TestScoreL);
CompCompleteTime = meanStdProcess(CompCompleteTime);
CompSolveTime = meanStdProcess(CompSolveTime);
CompScore = meanStdProcess(CompScore);
CompScoreL = meanStdProcess(CompScoreL);
degreeConn = meanStdProcess(degreeConn);
N = meanStdProcess(N);
% if ~isempty(equiPerformanceTime)
%     equiPerformanceTime = meanStdProcess(equiPerformanceTime);
% end
if ~isempty(equiPerformanceTimeL)
    equiPerformanceTimeL = meanStdProcess(equiPerformanceTimeL);
end
equiTimePerformance = meanStdProcess(equiTimePerformance);

TestTotalTime = meanStdProcess(TestTotalTime);
CompTotalTime = meanStdProcess(CompTotalTime);


%% plot
figure(1)
clf
mesh(N(:,:,1),degreeConn(:,:,1),(CompTotalTime(:,:,1)-TestTotalTime(:,:,1))./TestTotalTime(:,:,1)*100,'FaceColor',[1 0 0])
% hold on
% mesh(N(:,:,1),degreeConn(:,:,1),CompTotalTime(:,:,1))
hold on
mesh(N(:,:,1),degreeConn(:,:,1),zeros(size(N(:,:,1))))
xlabel('nodeNum')
ylabel('connectivity')
title('calc time')
view(-45,45)

figure(2)
clf
mesh(N(:,:,1),degreeConn(:,:,1),(CompScore(:,:,1)-TestScore(:,:,1))./TestScore(:,:,1)*100,'FaceColor',[1 0 0])
% hold on
% mesh(N(:,:,1),degreeConn(:,:,1),CompScore(:,:,1))
hold on
mesh(N(:,:,1),degreeConn(:,:,1),zeros(size(N(:,:,1))))
xlabel('nodeNum')
ylabel('connectivity')
title('test score')
view(-45,45)

figure(3)
clf
mesh(N(:,:,1),degreeConn(:,:,1),(CompScoreL(:,:,1)-TestScoreL(:,:,1))./TestScoreL(:,:,1)*100,'FaceColor',[1 0 0])
% hold on
% mesh(N(:,:,1),degreeConn(:,:,1),CompScoreL(:,:,1))
hold on
mesh(N(:,:,1),degreeConn(:,:,1),zeros(size(N(:,:,1))))
xlabel('nodeNum')
ylabel('connectivity')
title('test score L')
view(-45,45)

figure(4)
clf
mesh(N(:,:,1),degreeConn(:,:,1),equiPerformanceTimeL(:,:,1),'FaceColor',[1 0 0])
xlabel('nodeNum')
ylabel('connectivity')
title('equiPerfTimeL')
view(-45,45)

figure(5)
clf
mesh(N(:,:,1),degreeConn(:,:,1),equiPerformanceTime(:,:,1),'FaceColor',[1 0 0])
xlabel('nodeNum')
ylabel('connectivity')
title('equiPerfTimeM')
view(-45,45)

figure(6)
clf
mesh(N(:,:,1),degreeConn(:,:,1),(equiTimePerformance(:,:,1)-TestScoreL(:,:,1))./TestScoreL(:,:,1)*100,'FaceColor',[1 0 0])
% mesh(N(:,:,1),degreeConn(:,:,1),(equiTimePerformance(:,:,1)),'FaceColor',[1 0 0])
% hold on
% mesh(N(:,:,1),degreeConn(:,:,1),CompScoreL(:,:,1))
hold on
mesh(N(:,:,1),degreeConn(:,:,1),zeros(size(N(:,:,1))))
xlabel('nodeNum')
ylabel('connectivity')
title('equiTimePerf')
view(-45,45)

% end


function output = meanStdProcess(input)
m = mean(input,3);
s = std(input,[],3);
output(:,:,1) = m;
output(:,:,2) = s;
end