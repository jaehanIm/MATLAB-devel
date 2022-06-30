addpath('./data/')

%% Load data
dat = load('ML6.mat'); % ML4 is for paper
dat = dat.MCData;
dat2 = load('ML6_missing.mat');
dat2 = dat2.MCData;
dat3 = load('ML6_missing2.mat');
dat3 = dat3.MCData;
dat4 = load('ML6_missing3.mat');
dat4 = dat4.MCData;
% dat(end+1,:,:) = dat2;
% dat(end+1,:,:) = dat3;
dat(4,:,:) = [];
dat(end+1,:,:) = dat4;

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
            if sizeFactor == 8 & connFactor == 1 & i == 9
                dat(sizeFactor,connFactor,i) = dat(sizeFactor,connFactor,i-1);
            end
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

%%
degreeRange = [0,0.05,0.1,0.2,0.35,0.5];
filterDatTime = [];
filterDatScore = [];
for k = 1:size(degreeRange,2)-1
    temp1 = []; temp2 = [];
    for i = 1:sizeFactor*connFactor
        if degreeConn(i) > degreeRange(k) && degreeConn(i) < degreeRange(k+1)
            temp1 = vertcat(temp1,[N(i),scoreRatio(i)]);
            temp2 = vertcat(temp2,[N(i),timeRatio(i)]);
        end
    end
    filterDatScore{k} = temp1;
    filterDatTime{k} = temp2;
end
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
% figure(11)
% clf
% for i = 1:size(degreeConn,1)
% semilogx(degreeConn(i,:,1),CompTotalTime(i,:,1),'*--')
% hold on
% semilogx(degreeConn(i,:,1),TestTotalTime(i,:,1),'o-')
% end
% grid on
% xlabel('connectivity')
% ylabel('computation time [s]')
% title('Computation Time')

figure(112)
clf
for i = 1:size(degreeConn,1)
% h = errorbar(degreeConn(i,:,1),timeRatio(i,:,1),timeRatio(i,:,2),'-*');
semilogy(degreeConn(i,:,1),timeRatio(i,:,1),'-*');
text(degreeConn(i,end,1),timeRatio(i,end,1),num2str(Nlist(i)))
hold on
end
grid on
xlabel('connectivity')
ylabel('computation time [%]')
title('Computation Time Ratio')
% set(get(h,'Parent'), 'XScale','log')

% figure(113)
% clf
% for i = 1:size(degreeConn,1)
% plot(degreeConn(i,:,1),CompScoreL(i,:,1),'*--')
% hold on
% plot(degreeConn(i,:,1),TestScoreL(i,:,1),'o-')
% end
% grid on
% xlabel('connectivity')
% ylabel('score')
% title('Score')
% 
% figure(22)
% clf
% for i = 1:size(degreeConn,1)
% h = errorbar(degreeConn(i,:,1),scoreRatio(i,:,1),scoreRatio(i,:,2),'*');
% text(degreeConn(i,end,1),scoreRatio(i,end,1),num2str(Nlist(i)))
% hold on
% end
% grid on
% xlabel('connectivity')
% ylabel('score ratio [%]')
% title('score ratio ')
% % set(get(h,'Parent'), 'YScale','log')

% figure(33)
% clf
% temp = N(:,:,1);
% temp2 = CompTotalTime(:,:,1);
% temp3 = TestTotalTime(:,:,1);
% loglog(temp(:),temp2(:),'*')
% hold on
% loglog(temp(:),temp3(:),'*')
% grid on
% xlabel('Node Number')
% ylabel('Computation Time')
% title('Computation Time Comparison')
% 
% figure(332)
% clf
% temp = N(:,:,1);
% temp2 = CompScoreL(:,:,1);
% temp3 = TestScoreL(:,:,1);
% semilogx(temp(:),temp2(:),'*')
% hold on
% semilogx(temp(:),temp3(:),'*')
% grid on
% xlabel('Node Number')
% ylabel('Score')
% title('Score Comparison')

figure(333)
clf
temp = N(:,:,1);
temp2 = CompScoreL(:,:,1)./TestScoreL(:,:,1);
semilogx(temp(:),temp2(:),'*')
grid on
xlabel('Node Number')
ylabel('Score')
title('Score ratio')

figure(334)
clf
temp = N(:,:,1);
temp2 = CompTotalTime(:,:,1)./TestTotalTime(:,:,1);
loglog(temp(:),temp2(:),'*')
grid on
xlabel('Node Number')
ylabel('comp time ratio')
title('Computation time ratio')

% figure(441)
% clf
% temp = degreeConn(:,:,1);
% temp2 = CompTotalTime(:,:,1);
% temp3 = TestTotalTime(:,:,1);
% plot(temp(:),temp2(:),'*')
% hold on
% plot(temp(:),temp3(:),'*')
% grid on
% xlabel('connectivity')
% ylabel('Computation Time')
% title('Computation Time Comparison')
% 
% figure(442)
% clf
% temp = degreeConn(:,:,1);
% temp2 = CompScoreL(:,:,1);
% temp3 = TestScoreL(:,:,1);
% plot(temp(:),temp2(:),'*')
% hold on
% plot(temp(:),temp3(:),'*')
% grid on
% xlabel('connectivity')
% ylabel('Score')
% title('Score Comparison')

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

% figure(551)
% clf
% for i = 1:size(degreeRange,2)-1
%     loglog(filterDatTime{i}(:,1),filterDatTime{i}(:,2),'*')
%     hold on
% end
% grid on
% xlabel('Node num')
% ylabel('time ratio')
% legend('1','2','3','4','5','6')
% title('time ratio')
% 
figure(552)
clf
for i = 1:size(degreeRange,2)-1
    semilogx(filterDatScore{i}(:,1),filterDatScore{i}(:,2),'x','LineWidth',5)
    hold on
end
grid on
xlabel('Node num')
ylabel('score ratio')
legend('Low','Low-moderate','Moderate','High-moderate','High')
title('score ratio')


figure(1)
clf
temp = degreeConn(:,:,1);
temp1 = N(:,:,1);
temp2 = CompCompleteTime(:,:,1);
temp3 = TestCompleteTime(:,:,1);
temp4 = TestClusteringTime(:,:,1);
% semilogy(temp(:),temp2(:)./(temp3(:)+temp4(:)),'*')
for i = 1:size(degreeConn,1)
    semilogy(temp(i,:,1),temp2(i,:,1)./(temp3(i,:,1)+temp4(i,:,1)),'*-')
    text(temp(i,:,1),temp2(i,:,1)./(temp3(i,:,1)+temp4(i,:,1)),num2str(Nlist(i)))
    hold on
end
grid on
title('Complefication time ratio')
xlabel('connectivity')
figure(2)
clf
loglog(temp1(:),temp2(:)./(temp3(:)+temp4(:)),'*')
grid on
title('Complefication time ratio')
xlabel('Node number')
% figure(3)
% temp5 = CompSolveTime(:,:,1);
% temp6 = TestSolveTime(:,:,1);
% clf
% plot(temp1(:),temp5(:)./temp6(:),'*')
% grid on
% title('Solve time ratio')
% xlabel('node number')
figure(4)
temp5 = CompSolveTime(:,:,1);
temp6 = TestSolveTime(:,:,1);
clf
semilogx(temp1(:),temp5(:)./temp6(:),'*')
grid on
title('Solve time ratio')
xlabel('node Number')

figure(5)
clf
for i = 1:size(degreeConn,1)
color = rand(1,3);
h = plot(degreeConn(i,:,1),CompSolveTime(i,:,1)./TestSolveTime(i,:,1),'*-','Color',color);
text(degreeConn(i,end,1),CompSolveTime(i,end,1)./TestSolveTime(i,end,1),num2str(Nlist(i)))
hold on
end
grid on
xlabel('connectivity')
ylabel('solve time ratio')
title('solve time ratio')

figure(6)
clf
temp = N(:,:,1);
temp1 = CompSolveTime(:,:,1)./CompTotalTime(:,:,1)*100;
temp1_2 = CompCompleteTime(:,:,1)./CompTotalTime(:,:,1)*100;
temp2 = TestSolveTime(:,:,1)./TestTotalTime(:,:,1)*100;
temp2_2 = TestCompleteTime(:,:,1)./TestTotalTime(:,:,1)*100;
temp2_3 = TestClusteringTime(:,:,1)./TestTotalTime(:,:,1)*100;
for i = 1:size(degreeConn,1)
    semilogx(temp(:,i,1),temp1_2(:,i,1),'rx')
    hold on
    semilogx(temp(:,i,1),temp2_2(:,i,1)+temp2_3(:,i,1),'bo')
end
% for i = 1:size(degreeConn,1)
%     semilogx(temp(:,i,1),temp1(:,i,1),'rx')
%     hold on
%     semilogx(temp(:,i,1),temp2(:,i,1),'bo')
% end
grid on
xlabel('Node number')
ylabel('compostion ratio [%]')
title('The ratio of time spent on complete graph construction')
legend('Pure ACS', 'Proposed Algorithm')
%% function

function output = meanStdProcess(input)
m = mean(input,3,'omitnan');
s = std(input,[],3,'omitnan')/sqrt(size(input,3));
output(:,:,1) = m;
output(:,:,2) = s;
end
