% function MC_analysis()
addpath('./data/')

%% Load data
dat = load('set3.mat');
dat = dat.MCData;
dat(end,:,:) = [];
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

iterN = size(dat,3);

sizeList = 3:-0.5:0.5;
connList = 5:5:50;
[X,Y] = meshgrid(connList,sizeList);

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

TestTotalTime = meanStdProcess(TestTotalTime);
CompTotalTime = meanStdProcess(CompTotalTime);


%% plot
figure(1)
clf
mesh(N(:,:,1),degreeConn(:,:,1),TestTotalTime(:,:,1),'FaceColor',[1 0 0])
hold on
mesh(N(:,:,1),degreeConn(:,:,1),CompTotalTime(:,:,1))
xlabel('nodeNum')
ylabel('connectivity')

figure(2)
clf
mesh(N(:,:,1),degreeConn(:,:,1),TestScore(:,:,1),'FaceColor',[1 0 0])
hold on
mesh(N(:,:,1),degreeConn(:,:,1),CompScore(:,:,1))
xlabel('nodeNum')
ylabel('connectivity')

figure(3)
clf
mesh(N(:,:,1),degreeConn(:,:,1),TestScoreL(:,:,1),'FaceColor',[1 0 0])
hold on
mesh(N(:,:,1),degreeConn(:,:,1),CompScoreL(:,:,1))
xlabel('nodeNum')
ylabel('connectivity')




% end


function output = meanStdProcess(input)
m = mean(input,3);
s = std(input,[],3);
output(:,:,1) = m;
output(:,:,2) = s;
end