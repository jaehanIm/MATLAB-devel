addpath('./../')
addpath('..\ACO\')

%%%%%%%%%%% ACO for VRPTW version %%%%%%%%%%%

% parameter setting
fovFactor = 2.8;
mapheight = 3;
inpection_dist = 7;

distThres = 15;
vnum = 8;
antNo = 20;
stopThres = 200;
capacity = 395;

% generate map
mapGenerator_VRPTW
node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
node = vertcat([20,5,0],node);

% node = [0,0,0;10,10,0;20,0,0;-30,5,0;-25,-15,0]/2;

% load instance
% no. x. y. demand. init. due. servT.
% map = load('instance.mat');
% map = map.temp;
% node = map(:,2:3); node = horzcat(node,zeros(size(map,1),1));

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
A=A_orig;
toc

% mapData generation
mapGraph.n = N;
mapGraph.edges = C;
mapGraph.node = node;

% free timeWindow
timeWindow = zeros(N,2);
timeWindow(:,1) = 0;
timeWindow(:,2) = 10000;
servTime = zeros(N,1);

% time constraint construction
% [~,~,timeMax] = NNHeuristic_VRPTW(mapGraph);
% timeWindow = zeros(N,2);
% timeWindow(:,1) = 0;
% timeWindow(:,2) = timeMax;
% timeWindow = map(:,5:6);
% servTime = map(:,7)/2;
% nodeDemand = map(:,4);

% timeWindow = [0, 100; 2,30; 10,40; 60,80; 30,60];
% servTime = [10;10;10;10;10];


% solve problem
ACS_VRPTW;

% plot result
drawBestTour_forSoleVRPTW(colony, mapGraph, vnum);
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
title('VRPTW')

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

%% Animation

figure(5)
clf
% draw node
for i = 1:size(node,1)
    plot3(node(:,1),node(:,2),node(:,3),'.');
end
lineList = find(A(:));
for i = 1:N-1
    for j= i+1:N
        if A(i,j) ~= 0
            initIdx = i;
            termIdx = j;
            line([node(initIdx,1), node(termIdx,1)],[node(initIdx,2), node(termIdx,2)],[node(initIdx,3), node(termIdx,3)]);
        end
    end
end
hold on
plot3(node(1,1),node(1,2),node(1,3),'x','MarkerSize',5,'LineWidth',4)
grid on
axis equal
title('Tour animation')
% view(0, 90)

% animate
tour = colony.queen.tour;
tick = colony.queen.tickHistory;
tourLen = colony.queen.vehTourLen;
finished = zeros(vnum,1);

% vid = VideoWriter('animation.avi','MPEG-4');
% vid.FrameRate = fps;
% vid.Quality = 100;
vid = VideoWriter('animation_comparison.avi');
open(vid);

simStep = 1;
T = max(tick,[],'all');
stepNum = ceil(T/simStep);
simT = linspace(0,T,stepNum);

for i = 1:vnum
    vv(i) = plot3(node(1,1),node(1,2),node(1,3),'o','LineWidth',3);
    hold on
end
hh(i) = line([0,0],[0,0],[0,0]);
view(0,90)

count = 0;
for t = simT
    count = count + 1;
    for v = 1:vnum
        % calc status
        routeIdx = getRouteStep(tick(v,:),t);
        if routeIdx ~= -1
            initNode = tour(v,routeIdx); termNode = tour(v,routeIdx+1);
            if A(initNode, termNode) ~= 0 && termNode ~= 1
                regionDur = tick(v,routeIdx+1) - tick(v,routeIdx);
                weight = (t-tick(v,routeIdx))/regionDur;
                if t-tick(v,routeIdx) > regionDur
                    weight = 1;
                end
                newPos = (1-weight) .* node(initNode,:) + weight .* node(termNode,:);
                set(vv(v), 'XData', newPos(1), 'YData', newPos(2), 'ZData', newPos(3));
            elseif termNode == 1
                finished(v) = true;
                set(vv(v), 'XData', node(tour(v,tourLen(v)),1), 'YData', node(tour(v,tourLen(v)),2), 'ZData', node(tour(v,tourLen(v)),3));
            else
                bypassRoute = implicitRoute{initNode, termNode};
                occupancyInfo = occupancy{v,routeIdx+1}(:,3:4);
                occupancyInit = occupancyInfo(:,1);
                subRouteIdx = max(find((occupancyInit - t) < 0));
                if ~isempty(subRouteIdx) && min(occupancyInit) <= t
                    subRouteInfo = occupancyInfo(subRouteIdx,:);
                    regionDur = subRouteInfo(2) - subRouteInfo(1);
                    weight = (t - subRouteInfo(1))/regionDur;
                    if (t - subRouteInfo(1)) > regionDur
                        weight = 1;
                    end
                    newPos = (1-weight) .* node(bypassRoute(subRouteIdx),:) + weight .* node(bypassRoute(subRouteIdx+1),:);
                    set(vv(v), 'XData', newPos(1), 'YData', newPos(2), 'ZData', newPos(3));
                end
            end
        elseif finished(v) == true
        else
            set(vv(v), 'XData', node(1,1), 'YData', node(1,2), 'ZData', node(1,3));
        end
    end
%     for i = 1:N
%         for j = 1:N
%             reservationLen = reservation{i,j}.num;
%             if reservationLen ~= 0
%                 for k = 1:reservationLen
%                     timeRange = reservation{i,j}.info{k}(1:2);
%                     if (timeRange(1)-t) * (timeRange(2)-t) < 0
%                     end
%                 end
%             end
%         end
%     end
    drawnow;
    frame = getframe(gcf);
    writeVideo(vid,frame);
    pause(0.01)
end
close(vid);
