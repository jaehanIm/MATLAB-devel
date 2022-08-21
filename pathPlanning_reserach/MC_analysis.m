addpath('./data/')

%% Load data
dat = load('MLFF.mat'); % ML4 is for paper
dat = dat.MCData;
% dat = temp;
dat2 = load('MLF_missing.mat');
dat2 = dat2.MCData;
dat3 = load('MLF_missing2.mat');
dat3 = dat3.MCData;
dat4 = load('MLF_missing3.mat');
dat4 = dat4.MCData;
% dat(end+1,:,:) = dat2;
% dat(end+1,:,:) = dat3;
% dat(4,:,:) = [];
% dat(end+1,:,:) = dat4;

%% Additional data - MY BIG MISTAKE
temp = load('MLF_supp.mat');
temp = temp.MCData;
temp2 = load('MLF_supp_missing1.mat');
temp2 = temp2.MCData;
temp3 = load('MLF_supp_missing2.mat');
temp3 = temp3.MCData;
temp4 = load('MLF_supp_missing3.mat');
temp4 = temp4.MCData;
% temp(end+1,:) = temp2;
% temp(end+1,:) = temp3;
% temp(4,:) = [];
% temp(end+1,:) = temp4;

% dat(:,5:end,:) = [];

for i = 1:size(dat,1)
    for j = 1:size(dat,2)
        for k = 1:20
            dat{i,j,k}.TestIntraCompleteTime = temp{i,j}.TestIntraCompleteTime;
        end
    end
end

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
%             if sizeFactor == 8 & connFactor == 1 & i == 19
%                 dat(sizeFactor,connFactor,i) = dat(sizeFactor,connFactor,i-1);
%             end
            if sizeFactor == 9 & connFactor == 1 & i == 19
                dat(sizeFactor,connFactor,i) = dat(sizeFactor,connFactor,i-1);
            end
            TestClusteringTime(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.TestClusteringTime;
            TestCompleteTime(sizeFactor,connFactor,i) = dat{sizeFactor,connFactor,i}.TestIntraCompleteTime + dat{sizeFactor,connFactor,i}.TestInterCompleteTime;
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

% TestClusteringTime = zeros(size(TestClusteringTime));
TestTotalTime = TestClusteringTime + TestCompleteTime + TestSolveTime;
CompTotalTime = CompCompleteTime + CompSolveTime;
% timeRatio = (CompTotalTime - TestTotalTime)./TestTotalTime*100;
timeRatio = CompTotalTime./TestTotalTime;
completeTimeRatio = CompCompleteTime./TestCompleteTime;
solveTimeRatio = CompSolveTime./TestSolveTime;
preTimeRatio = CompCompleteTime./(TestCompleteTime + TestClusteringTime);
% scoreRatio = (CompScoreL - TestScoreL)./TestScoreL*100;
scoreRatio = CompScoreL./TestScoreL;
TestTotalTimeRaw = TestTotalTime;
CompTotalTimeRaw = CompTotalTime;
timeRatioRaw = timeRatio;
scoreRatioRaw = scoreRatio;
CompScoreLRaw = CompScoreL;
TestScoreLRaw = TestScoreL;
completeTimeRatioRaw = completeTimeRatio;
solveTimeRatioRaw = solveTimeRatio;
preTimeRatioRaw = preTimeRatio;

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
completeTimeRatio = meanStdProcess(completeTimeRatio);
solveTimeRatio = meanStdProcess(solveTimeRatio);
preTimeRatio = meanStdProcess(preTimeRatio);

%%
% degreeRange = [0,0.05,0.1,0.2,0.35,0.5];
degreeRange = [0,0.1,0.3,0.5];
% degreeRange = [0,0.05,0.1];
filterDatTime = []; filterDatScore = []; filterDatCompleteTime = []; filterDatSolveTime = []; filterDatPreTime = [];
for k = 1:size(degreeRange,2)-1
    temp1 = [];  temp2 = []; temp3 = []; temp4 = []; temp5 = [];
    for i = 1:sizeFactor*connFactor
        if degreeConn(i) > degreeRange(k) && degreeConn(i) < degreeRange(k+1)
            temp1 = vertcat(temp1,[N(i),scoreRatio(i)]);
            temp2 = vertcat(temp2,[N(i),timeRatio(i)]);
            temp3 = vertcat(temp3,[N(i),completeTimeRatio(i)]);
            temp4 = vertcat(temp4,[N(i),solveTimeRatio(i)]);
            temp5 = vertcat(temp5,[N(i),preTimeRatio(i)]);
        end
    end
    filterDatScore{k} = temp1;
    filterDatTime{k} = temp2;
    filterDatCompleteTime{k} = temp3;
    filterDatSolveTime{k} = temp4;
    filterDatPreTime{k} = temp5;
end

NRange = [0,100,1000,4000];
NfilterDatTime = []; NfilterDatScore = []; NfilterDatCompleteTime = []; NfilterDatSolveTime = []; NfilterDatPreTime = [];
for k = 1:size(NRange,2)-1
    temp1 = []; temp2 = []; temp3 = []; temp4 = []; temp5 = [];
    for i = 1:sizeFactor*connFactor
        if N(i) > NRange(k) && N(i) < NRange(k+1)
            temp1 = vertcat(temp1,[degreeConn(i),scoreRatio(i)]);
            temp2 = vertcat(temp2,[degreeConn(i),timeRatio(i)]);
            temp3 = vertcat(temp3,[degreeConn(i),completeTimeRatio(i)]);
            temp4 = vertcat(temp4,[degreeConn(i),solveTimeRatio(i)]);
            temp5 = vertcat(temp5,[degreeConn(i),preTimeRatio(i)]);
        end
    end
    NfilterDatScore{k} = temp1;
    NfilterDatTime{k} = temp2;
    NfilterDatCompleteTime{k} = temp3;
    NfilterDatSolveTime{k} = temp4;
    NfilterDatPreTime{k} = temp5;
end

%% mono plot
% legendList = [];
% for i = 1:10
%     legendList = vertcat(legendList," " + [num2str(i)]);
% end
Nlist = N(:,1,1);
shapeList = ["+","x","o","d","h"];
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
for i = 1:size(NRange,2)-1
    plot(NfilterDatScore{i}(:,1),NfilterDatTime{i}(:,2),shapeList(i),'LineWidth',1,'MarkerSize',7)
    hold on
end
for i = 1:size(degreeConn,1)
% h = errorbar(degreeConn(i,:,1),timeRatio(i,:,1),timeRatio(i,:,2),'k-*');
plot(degreeConn(i,:,1),timeRatio(i,:,1),':','Color',[0.6 0.6 0.6]);
% text(degreeConn(i,end,1)+0.01,timeRatio(i,end,1),num2str(Nlist(i)))
hold on
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
grid on
xlabel('$\mathcal{K}$','Interpreter','latex')
ylabel('CT ratio [log scale]')
title('CT Ratio')
% set(get(h,'Parent'), 'XScale','log')
legend('N <= 100','N <= 1000','N > 1000','Location','best')
set(gca,'children',flipud(get(gca,'children')))

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

% figure(333)
% clf
% temp = N(:,:,1);
% temp2 = CompScoreL(:,:,1)./TestScoreL(:,:,1);
% semilogx(temp(:),temp2(:),'*')
% grid on
% xlabel('Node Number')
% ylabel('Score')
% title('Score ratio')

figure(334)
clf
temp = N(:,:,1);
temp2 = CompTotalTime(:,:,1)./TestTotalTime(:,:,1);
% loglog(temp(:),temp2(:),'*')
for i = 1:size(degreeRange,2)-1
    semilogx(filterDatScore{i}(:,1),filterDatTime{i}(:,2),shapeList(i),'LineWidth',1,'MarkerSize',7)
    hold on
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
grid on
xlabel('$\textbf{N}$ [log scale]','Interpreter','latex')
ylabel('CT ratio [log scale]')
title('CT ratio')
legend('Low','Low-moderate','Moderate','High-moderate','High','Location','northwest')

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
% temp = degreeConn(:,:,1);
% temp2 = CompScoreL(:,:,1)./TestScoreL(:,:,1);
% plot(temp(:),temp2(:),'*')
for i = 1:size(NRange,2)-1
    plot(NfilterDatScore{i}(:,1),NfilterDatScore{i}(:,2),shapeList(i),'LineWidth',1,'MarkerSize',7)
    hold on
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
grid on
xlabel('$\mathcal{K}$','Interpreter','latex')
ylabel('SQ ratio')
title('SQ ratio')
legend('N <= 100','N <= 1000','N > 1000')

% figure(444)
% clf
% temp = degreeConn(:,:,1);
% temp2 = CompTotalTime(:,:,1)./TestTotalTime(:,:,1);
% plot(temp(:),temp2(:),'*')
% grid on
% xlabel('connectivity')
% ylabel('comp time ratio')
% title('Computation time ratio')

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
    semilogx(filterDatScore{i}(:,1),filterDatScore{i}(:,2),shapeList(i),'MarkerSize',7,'LineWidth',1)
    hold on
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
grid on
xlabel('$\textbf{N}$ [log scale]','Interpreter','latex')
ylabel('SQ ratio')
legend('Low','Low-moderate','Moderate','High-moderate','High')
title('SQ ratio')


figure(1)
clf
temp = degreeConn(:,:,1);
temp1 = N(:,:,1);
temp2 = CompCompleteTime(:,:,1);
temp3 = TestCompleteTime(:,:,1);
temp4 = TestClusteringTime(:,:,1);
for i = 1:size(NRange,2)-1
%     semilogy(NfilterDatCompleteTime{i}(:,1),NfilterDatCompleteTime{i}(:,2),shapeList(i),'LineWidth',1,'MarkerSize',7)
    semilogy(NfilterDatPreTime{i}(:,1),NfilterDatPreTime{i}(:,2),shapeList(i),'LineWidth',1,'MarkerSize',7)
    hold on
end
for i = 1:size(degreeConn,1)
    semilogy(temp(i,:,1),temp2(i,:,1)./(temp3(i,:,1)+temp4(i,:,1)),':','Color',[0.6 0.6 0.6]);
%     semilogy(temp(i,:,1),temp2(i,:,1)./temp3(i,:,1),':','Color',[0.6 0.6 0.6]);
    hold on
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
grid on
xlabel('$\mathcal{K}$','Interpreter','latex')
ylabel('CT_{pre} ratio [log scale]')
title('CT_{pre} ratio')
% set(get(h,'Parent'), 'XScale','log') 
legend('N <= 100','N <= 1000','N > 1000','Location','best')
set(gca,'children',flipud(get(gca,'children')))


figure(2)
clf
% loglog(temp1(:),temp2(:)./(temp3(:)+temp4(:)),'*')
for i = 1:size(degreeRange,2)-1
    loglog(filterDatPreTime{i}(:,1),filterDatPreTime{i}(:,2),shapeList(i),'LineWidth',1,'MarkerSize',7);
    hold on
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
grid on
title('CT_{pre} ratio')
xlabel('$\textbf{N}$ [log scale]','Interpreter','latex')
ylabel('CT_{pre} ratio [log scale]')
legend('Low','Low-moderate','Moderate','High-moderate','High','Location','northwest')

figure(3)
clf
% plot(temp(:),temp5(:)./temp6(:),'*')
for i = 1:size(NRange,2)-1
    plot(NfilterDatSolveTime{i}(:,1),NfilterDatSolveTime{i}(:,2),shapeList(i),'LineWidth',1,'MarkerSize',7)
    hold on
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
grid on
title('CT_{s} ratio')
xlabel('$\mathcal{K}$','Interpreter','latex')
ylabel('CT_{s} ratio')
legend('N <= 100','N <= 1000','N > 1000')

figure(4)
temp5 = CompSolveTime(:,:,1);
temp6 = TestSolveTime(:,:,1);
clf
% semilogx(temp1(:),temp5(:)./temp6(:),'*')
for i = 1:size(degreeRange,2)-1
    semilogx(filterDatSolveTime{i}(:,1),filterDatSolveTime{i}(:,2),shapeList(i),'LineWidth',1,'MarkerSize',7)
    hold on
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
grid on
title('CT_{s} ratio')
ylabel('CT_{s} ratio')
xlabel('$\textbf{N}$ [log scale]','Interpreter','latex')
legend('Low','Low-moderate','Moderate','High-moderate','High','Location','northwest')


% figure(5)
% clf
% for i = 1:size(degreeConn,1)
% % color = rand(1,3);
% h = plot(degreeConn(i,:,1),CompSolveTime(i,:,1)./TestSolveTime(i,:,1),'*-');
% text(degreeConn(i,end,1)+0.01,CompSolveTime(i,end,1)./TestSolveTime(i,end,1),num2str(Nlist(i)))
% hold on
% end
% grid on
% xlabel('connectivity')
% ylabel('solve time ratio')
% title('solve time ratio')

figure(6)
clf
temp = N(:,:,1);
temp1 = CompSolveTime(:,:,1)./CompTotalTime(:,:,1)*100;
temp1_2 = CompCompleteTime(:,:,1)./CompTotalTime(:,:,1)*100;
temp2 = TestSolveTime(:,:,1)./TestTotalTime(:,:,1)*100;
temp2_2 = TestCompleteTime(:,:,1)./TestTotalTime(:,:,1)*100;
temp2_3 = TestClusteringTime(:,:,1)./TestTotalTime(:,:,1)*100;
for i = 1:size(degreeConn,2)
    semilogx(temp(:,i,1),temp1_2(:,i,1),'rx')
    hold on
    semilogx(temp(:,i,1),temp2_2(:,i,1)+temp2_3(:,i,1),'bo')
end
lim = axis;
plot([lim(1) lim(2)],[1 1],'k:','LineWidth',2)
% for i = 1:size(degreeConn,1)
%     semilogx(temp(:,i,1),temp1(:,i,1),'rx')
%     hold on
%     semilogx(temp(:,i,1),temp2(:,i,1),'bo')
% end
grid on
xlabel('$\textbf{N}$ [log scale]','Interpreter','latex')
ylabel('proportion [%]')
% title('Ratio of time spent on complete graph construction')
% title('Ratio of CT_{cg}/CT')
title('Proportion of CT_{pre} to CT')
legend('Pure ACS', 'Proposed Algorithm','Location','northwest')
%% function

function output = meanStdProcess(input)
m = mean(input,3,'omitnan');
s = std(input,[],3,'omitnan')/sqrt(size(input,3));
output(:,:,1) = m;
output(:,:,2) = s;
end
