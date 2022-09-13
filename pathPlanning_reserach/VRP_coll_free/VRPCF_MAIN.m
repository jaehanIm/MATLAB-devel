clear all

addpath('./../')
addpath('..\ACO\')

%%%%%%%%%%% ACO for VRPCF version %%%%%%%%%%%

% parameter setting
fovFactor = 3;
mapheight = 3;
inpection_dist = 7;

distThres = 25;
vnum = 6;
antNo = 20;
stopThres = 50;
capacity = 395;

homePos = [20,20,6];

% generate map
mapGenerator_VRPCF
node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
node = vertcat(homePos,node);

N = size(node,1);
servTime = zeros(N,1);

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
C_orig = C;
tic
implicitRoute = cell(N,N);
for i = 1:N-1
    for j = 2:N
        if A(i,j) == 0
            [implRoute, implCost, ~] = shortestpath(graph(C_orig),i,j);
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
mapGraph.A = A;

% solve problem
ACS_VRPCF;

% plot
drawSchedule(2, N, colony, vnum);

%%

figure(4)
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
reservation = colony.queen.reservation;
occupancy = colony.queen.occupancy;

vid = VideoWriter('animation.avi');
open(vid);

simStep = 1;
T = max(tick,[],'all');
stepNum = ceil(T/simStep);
simT = linspace(0,T,stepNum);
mov(1:stepNum) = struct('cdata', [], 'colormap', []);

for i = 1:vnum
    vv(i) = plot3(node(1,1),node(1,2),node(1,3),'o','LineWidth',3);
    hold on
end
hh(i) = line([0,0],[0,0],[0,0]);

count = 0;
for t = simT
    count = count + 1;
    for v = 1:vnum
        % calc status
        routeIdx = getRouteStep(tick(v,:),t);
        if routeIdx ~= -1
            initNode = tour(v,routeIdx); termNode = tour(v,routeIdx+1);
            if A(initNode, termNode) ~= 0 || termNode == 1
                regionDur = tick(v,routeIdx+1) - tick(v,routeIdx);
                weight = (t-tick(v,routeIdx))/regionDur;
                newPos = (1-weight) .* node(initNode,:) + weight .* node(termNode,:);
                set(vv(v), 'XData', newPos(1), 'YData', newPos(2), 'ZData', newPos(3));
            else
                bypassRoute = implicitRoute{initNode, termNode};
                occupancyInfo = occupancy{v,routeIdx+1}(:,3:4);
                occupancyInit = occupancyInfo(:,1);
                subRouteIdx = max(find((occupancyInit - t)<0));
                if ~isempty(subRouteIdx)
                    subRouteInfo = occupancyInfo(subRouteIdx,:);
                    regionDur = subRouteInfo(2) - subRouteInfo(1);
                    weight = (t - subRouteInfo(1))/regionDur;
                    newPos = (1-weight) .* node(bypassRoute(subRouteIdx),:) + weight .* node(bypassRoute(subRouteIdx+1),:);
                    set(vv(v), 'XData', newPos(1), 'YData', newPos(2), 'ZData', newPos(3));                    
                end
            end
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

%%

figure(1)
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
