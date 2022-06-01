addpath('./data/')

%% Load data
dat = load('ML4.mat');
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
            if ~isempty(dat{sizeFactor,connFactor,i}.equiTimePerformance)
                equiTimePerformance(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.equiTimePerformance;
            else
                equiTimePerformance(sizeFactor,connFactor,i) = 0;
            end
        end
    end
end

TestTotalTime = TestClusteringTime + TestCompleteTime + TestSolveTime;
CompTotalTime = CompCompleteTime + CompSolveTime;
% timeRatio = (CompTotalTime - TestTotalTime)./TestTotalTime*100;
timeRatio = CompTotalTime./TestTotalTime;
% scoreRatio = (CompScoreL - TestScoreL)./TestScoreL*100;
scoreRatio = CompScoreL./TestScoreL;
TestTotalTimeRaw = TestTotalTime;
CompTotalTimeRaw = CompTotalTime;
timeRatioRaw = timeRatio;
scoreRatioRaw = scoreRatio;
CompScoreLRaw = CompScoreL;
TestScoreLRaw = TestScoreL;

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

if ~isempty(equiPerformanceTimeL)
    equiPerformanceTimeL = meanStdProcess(equiPerformanceTimeL);
end
equiTimePerformance = meanStdProcess(equiTimePerformance);

TestTotalTime = meanStdProcess(TestTotalTime);
CompTotalTime = meanStdProcess(CompTotalTime);
timeRatio = meanStdProcess(timeRatio);
scoreRatio = meanStdProcess(scoreRatio);


%% plot
% figure(1)
% clf
% mesh(N(:,:,1),degreeConn(:,:,1),(CompTotalTime(:,:,1)-TestTotalTime(:,:,1))./TestTotalTime(:,:,1)*100,'FaceColor',[1 0 0])
% % hold on
% % mesh(N(:,:,1),degreeConn(:,:,1),CompTotalTime(:,:,1))
% hold on
% mesh(N(:,:,1),degreeConn(:,:,1),zeros(size(N(:,:,1))))
% xlabel('nodeNum')
% ylabel('connectivity')
% title('calc time - red high better')
% view(-45,45)
% 
% figure(2)
% clf
% mesh(N(:,:,1),degreeConn(:,:,1),(CompScore(:,:,1)-TestScore(:,:,1))./TestScore(:,:,1)*100,'FaceColor',[1 0 0])
% % hold on
% % mesh(N(:,:,1),degreeConn(:,:,1),CompScore(:,:,1))
% hold on
% mesh(N(:,:,1),degreeConn(:,:,1),zeros(size(N(:,:,1))))
% xlabel('nodeNum')
% ylabel('connectivity')
% title('test score')
% view(-45,45)
% 
% figure(3)
% clf
% mesh(N(:,:,1),degreeConn(:,:,1),(CompScoreL(:,:,1)-TestScoreL(:,:,1))./TestScoreL(:,:,1)*100,'FaceColor',[1 0 0])
% % hold on
% % mesh(N(:,:,1),degreeConn(:,:,1),CompScoreL(:,:,1))
% hold on
% mesh(N(:,:,1),degreeConn(:,:,1),zeros(size(N(:,:,1))))
% xlabel('nodeNum')
% ylabel('connectivity')
% title('test score L - red high better')
% view(-45,45)
% 
% figure(4)
% clf
% mesh(N(:,:,1),degreeConn(:,:,1),equiPerformanceTimeL(:,:,1),'FaceColor',[1 0 0])
% xlabel('nodeNum')
% ylabel('connectivity')
% title('equiPerfTimeL - red high better')
% view(-45,45)
% 
% figure(5)
% clf
% mesh(N(:,:,1),degreeConn(:,:,1),equiPerformanceTime(:,:,1),'FaceColor',[1 0 0])
% xlabel('nodeNum')
% ylabel('connectivity')
% title('equiPerfTime')
% view(-45,45)
% 
% figure(6)
% clf
% mesh(N(:,:,1),degreeConn(:,:,1),(equiTimePerformance(:,:,1)-TestScoreL(:,:,1))./TestScoreL(:,:,1)*100,'FaceColor',[1 0 0])
% % mesh(N(:,:,1),degreeConn(:,:,1),(equiTimePerformance(:,:,1)),'FaceColor',[1 0 0])
% % hold on
% % mesh(N(:,:,1),degreeConn(:,:,1),CompScoreL(:,:,1))
% hold on
% mesh(N(:,:,1),degreeConn(:,:,1),zeros(size(N(:,:,1))))
% xlabel('nodeNum')
% ylabel('connectivity')
% title('equiTimePerf')
% view(-45,45)

%% mono plot
% legendList = [];
% for i = 1:10
%     legendList = vertcat(legendList," " + [num2str(i)]);
% end
Nlist = N(:,1,1);
figure(11)
clf
for i = 1:size(degreeConn,1)
semilogx(degreeConn(i,:,1),CompTotalTime(i,:,1),'*--')
hold on
semilogx(degreeConn(i,:,1),TestTotalTime(i,:,1),'o-')
end
grid on
xlabel('connectivity')
ylabel('computation time [s]')
title('Computation Time')

figure(112)
clf
for i = 1:size(degreeConn,1)
h = errorbar(degreeConn(i,:,1),timeRatio(i,:,1),timeRatio(i,:,2),'*--');
text(degreeConn(i,end,1),timeRatio(i,end,1),num2str(Nlist(i)))
hold on
end
grid on
xlabel('connectivity')
ylabel('computation time [%]')
title('Computation Time Ratio')
% set(get(h,'Parent'), 'XScale','log')

figure(113)
clf
for i = 1:size(degreeConn,1)
plot(degreeConn(i,:,1),CompScoreL(i,:,1),'*--')
hold on
plot(degreeConn(i,:,1),TestScoreL(i,:,1),'o-')
end
grid on
xlabel('connectivity')
ylabel('score')
title('Score')

figure(22)
clf
for i = 1:size(degreeConn,1)
h = errorbar(degreeConn(i,:,1),scoreRatio(i,:,1),scoreRatio(i,:,2),'*--');
text(degreeConn(i,end,1),scoreRatio(i,end,1),num2str(Nlist(i)))
hold on
end
grid on
xlabel('connectivity')
ylabel('score ratio [%]')
title('score ratio ')
% set(get(h,'Parent'), 'YScale','log')

figure(33)
clf
temp = N(:,:,1);
temp2 = CompTotalTime(:,:,1);
temp3 = TestTotalTime(:,:,1);
plot(temp(:),temp2(:),'*')
hold on
plot(temp(:),temp3(:),'*')
grid on
xlabel('Node Number')
ylabel('Computation Time')
title('Computation Time Comparison')

figure(332)
clf
temp = N(:,:,1);
temp2 = CompScoreL(:,:,1);
temp3 = TestScoreL(:,:,1);
plot(temp(:),temp2(:),'*')
hold on
plot(temp(:),temp3(:),'*')
grid on
xlabel('Node Number')
ylabel('Score')
title('Score Comparison')

figure(333)
clf
temp = N(:,:,1);
temp2 = CompScoreL(:,:,1)./TestScoreL(:,:,1);
plot(temp(:),temp2(:),'*')
grid on
xlabel('Node Number')
ylabel('Score')
title('Score ratio')

figure(334)
clf
temp = N(:,:,1);
temp2 = CompTotalTime(:,:,1)./TestTotalTime(:,:,1);
plot(temp(:),temp2(:),'*')
grid on
xlabel('Node Number')
ylabel('comp time ratio')
title('Computation time ratio')

% figure(3331)
% clf
% temp = N(:,:,1);
% temp2 = CompTotalTimeRaw(:,:,1);
% temp3 = TestTotalTimeRaw(:,:,1);
% plot(temp(:),temp2(:),'*')
% hold on
% plot(temp(:),temp3(:),'*')
% grid on
% xlabel('Node Number')
% ylabel('Computation Time')
% title('Computation Time Comparison')
% 
% figure(3332)
% clf
% temp = N(:,:,1);
% temp2 = CompScoreLRaw(:,:,1);
% temp3 = TestScoreLRaw(:,:,1);
% plot(temp(:),temp2(:),'*')
% hold on
% plot(temp(:),temp3(:),'*')
% grid on
% xlabel('Node Number')
% ylabel('Score')
% title('Score Comparison')
% 
% figure(3333)
% clf
% temp = N(:,:,1);
% temp2 = CompScoreLRaw(:,:,1)./TestScoreLRaw(:,:,1);
% plot(temp(:),temp2(:),'*')
% grid on
% xlabel('Node Number')
% ylabel('Score')
% title('Score ratioh')

figure(441)
clf
temp = degreeConn(:,:,1);
temp2 = CompTotalTime(:,:,1);
temp3 = TestTotalTime(:,:,1);
plot(temp(:),temp2(:),'*')
hold on
plot(temp(:),temp3(:),'*')
grid on
xlabel('connectivity')
ylabel('Computation Time')
title('Computation Time Comparison')

figure(442)
clf
temp = degreeConn(:,:,1);
temp2 = CompScoreL(:,:,1);
temp3 = TestScoreL(:,:,1);
plot(temp(:),temp2(:),'*')
hold on
plot(temp(:),temp3(:),'*')
grid on
xlabel('connectivity')
ylabel('Score')
title('Score Comparison')

figure(443)
clf
temp = degreeConn(:,:,1);
temp2 = CompScoreL(:,:,1)./TestScoreL(:,:,1);
plot(temp(:),temp2(:),'*')
grid on
xlabel('connectivity')
ylabel('Score')
title('Score ratio')

figure(444)
clf
temp = degreeConn(:,:,1);
temp2 = CompTotalTime(:,:,1)./TestTotalTime(:,:,1);
plot(temp(:),temp2(:),'*')
grid on
xlabel('connectivity')
ylabel('comp time ratio')
title('Computation time ratio')


%% function

function output = meanStdProcess(input)
m = mean(input,3);
s = std(input,[],3)/sqrt(size(input,3));
output(:,:,1) = m;
output(:,:,2) = s;
end
