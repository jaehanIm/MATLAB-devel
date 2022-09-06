addpath('./../')
addpath('..\ACO\')

%%%%%%%%%%% ACO for VRPCF version %%%%%%%%%%%

% parameter setting
fovFactor = 4;
mapheight = 3;
inpection_dist = 7;

distThres = 50;
vnum = 3;
antNo = 20;
stopThres = 200;
capacity = 395;

% generate map
mapGenerator_VRPCF
node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
node = vertcat([20,5,0],node);

N = size(node,1);

% Graph construction
L = zeros(N,N);
for i = 1:N-1
    for j = i+1:N
        L(i,j) = norm(node(i,:)-node(j,:));
        L(j,i) = L(i,j);
    end
end

A = zeros(N,N);
C = zeros(N,N);
for i = 1:N-1
    for j = 2:N
        if L(i,j) < distThres
            A(i,j) = 1;
            A(j,i) = 1;
            C(i,j) = L(i,j);
            C(j,i) = C(i,j);
        end
        if i == j
            A(i,j) = 0;
            A(j,i) = 0;
        end
    end
end

[A,C]=graphSparseConnection(node,A,C,L);
A_orig = A;
tic
implicitRoute = cell(N,N);
for i = 1:N-1
    for j = 2:N
        if A(i,j) == 0
            [implRoute, implCost, ~] = shortestpath(graph(C),i,j);
            implicitRoute{i,j} = implRoute;
            implicitRoute{j,i} = fliplr(implRoute);
            A(i,j) = 1; A(j,i) = 1;
            C(i,j) = implCost; C(j,i) = implCost;
        end
    end
end
toc
A = A_orig;

% mapData generation
mapGraph.n = N;
mapGraph.edges = C;
mapGraph.node = node;

% solve problem
ACS_VRPCF;

% plot result
drawBestTour_forSoleVRPCF(colony, mapGraph, vnum);
figure(2)
clf
hold on
for i = 1:N
    plot([timeWindow(i,1),timeWindow(i,2)],[i, i],'LineWidth',5,'Color',[0.3 0.3 0.3]')
end
grid on
ylim([1 N])
xlabel('time[s]')
ylabel('job number')
title('VRPCF')

plotTickHistory = colony.queen.tickHistory;
plotTour = colony.queen.tour;
plotTickHistory(:,1) = [];
plotTour(:,1) = [];

for v = 1:vnum
    color = rand(3,1)/2+0.5;
    for i = 1:length(plotTickHistory(v,plotTickHistory(v,:)~=0))
        startIdx = plotTour(v,i);
        endIdx = plotTour(v,i+1);
        plot([plotTickHistory(v,i)-servTime(startIdx), plotTickHistory(v,i)],[startIdx startIdx],'b','LineWidth',3,'Color',color)
    end
end

totTourLen = 0;
tourLen = zeros(vnum,1);
for v = 1:vnum
    for i = 1:length(colony.queen.tour(v,colony.queen.tour(v,:)~=0))-1
        startIdx = colony.queen.tour(v,i);
        endIdx = colony.queen.tour(v,i+1);
        tourLen(v) = tourLen(v) + C(startIdx,endIdx);
    end
end
totTourLen = sum(tourLen)

colony.queen.violation
sum(colony.queen.violation)