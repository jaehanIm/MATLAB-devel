clear all

addpath('./../')
addpath('..\ACO\')
addpath('./../Model/')

%%%%%%%%%%% ACO for VRPCF version %%%%%%%%%%%

% parameter setting
fovFactor = 2.8;
mapheight = 3;
inpection_dist = 7;

distThres = 20;
vnum = 4;
antNo = 20;
stopThres = 100;
capacity = 395;
servTime = 1;

% homePos = [30,40,6];
homePos = [10, 80, 6];
% homePos = [0,0,0];

global wowCount;
wowCount = 0;

% generate map
mapGenerator_VRPCF
node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
node = vertcat(homePos,node);
simStep = 1;

% temp node generator
% node = [0,0;1,1;1,-1;2,0;3,0;4,1;4,-1;5,0;1,2;1,3.1;2,2.5;2,1.5];
% node = horzcat(node,zeros(size(node,1),1));
% node = vertcat([1,3.9,0],node);
% vnum = 5; 
% distThres = sqrt(2)+0.02;
% simStep = 0.03;

% random node generator
% N = 70;
% node = rand(N,3);
% node(:,1:2) = node(:,1:2) * 100;
% node(:,3) = node(:,3) * 10;
% node = vertcat(homePos,node);
% simStep = 1;
% vnum = 15;

% random bridge generator
% N = 40;
% node = rand(N,3);
% node(:,1:2) = node(:,1:2) * 100 - 50;
% node2 = rand(N,3);
% node2(:,1:2) = node2(:,1:2) * 100 + 50;
% node = vertcat(node,node2);
% node(:,3) = node(:,3) * 10;
% node = vertcat(homePos,node);
% simStep = 1;
% vnum = 9;

% stl
% stlAddr = 'generic.stl';
% inpection_dist = 1;
% node = loadStl(stlAddr,inpection_dist);
% rndF = rand(size(node,1),1);
% node(rndF>0.3,:,:) = [];
% 
% % node segmentation
% body = node;
% 
% tail = body(body(:,1) > 14.7421 & body(:,2) < 0.955 & body(:,2) > -0.955,:,:);
% body(body(:,1) > 14.7421 & body(:,2) < 0.955 & body(:,2) > -0.955,:,:) = [];
% factor = 0.1; localN = size(tail,1);
% tail = tail(rand(localN,1) < factor,:,:);
% 
% lEngine = body(body(:,2) < -1.9 & body(:,2) > -5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) < 0.22,:,:);
% body(body(:,2) < -1.9 & body(:,2) > -5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) < 0.22,:,:) = [];
% factor = 0.1; localN = size(lEngine,1);
% lEngine = lEngine(rand(localN,1) < factor,:,:);
% 
% rEngine = body(body(:,2) > 1.9 & body(:,2) < 5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) < 0.22,:,:);
% body(body(:,2) > 1.9 & body(:,2) < 5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) < 0.22,:,:) = [];
% factor = 0.1; localN = size(rEngine,1);
% rEngine = rEngine(rand(localN,1) < factor,:,:);
% 
% emphennage = body(body(:,1)>12 & body(:,2)<2 & body(:,2)>-2,:,:);
% body(body(:,1)>12 & body(:,2)<2 & body(:,2)>-2,:,:) = [];
% factor = 0.2; localN = size(emphennage,1);
% emphennage = emphennage(rand(localN,1) < factor,:,:);
% 
% lEngineNeg = body(body(:,2) < -1.9 & body(:,2) > -5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) > 0.22 & body(:,3) < 0.5,:,:);
% body(body(:,2) < -1.9 & body(:,2) > -5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) > 0.22 & body(:,3) < 0.3,:,:) = [];
% 
% rEngineNeg = body(body(:,2) > 1.9 & body(:,2) < 5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) > 0.22 & body(:,3) < 0.5,:,:);
% body(body(:,2) > 1.9 & body(:,2) < 5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) > 0.22 & body(:,3) < 0.3,:,:) = [];
% 
% body(body(:,2) > -3.59 & body(:,2) < -3.51 & body(:,1) > 1.5 & body(:,1) < 3.95 & body(:,3) > 0.22,:,:) = [];
% body(body(:,2) > 3.51 & body(:,2) < 3.59 & body(:,1) > 1.5 & body(:,1) < 3.95 & body(:,3) > 0.22,:,:) = [];
% 
% nose = body(body(:,1) < -2.9,:,:);
% body(body(:,1) < -2.9,:,:) = [];
% factor = 0.5; localN = size(nose,1);
% nose = nose(rand(localN,1) < factor,:,:);
% 
% node = [body;tail;lEngine;rEngine;nose;emphennage];
% 
% save('stlnode.mat','node');
% 
% figure(98)
% clf
% drawStl(stlAddr,98)
% hold on
% plot3(node(:,1),node(:,2),node(:,3),'.')
% 
% simStep = 0.1;
% vnum = 9;
% distThres = 0.5;
% stopThres = 30;


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

graphDensity = sum(A,'all')/2/nchoosek(N,2)*100

% solve problem
ACS_VRPCF;

% plot
drawSchedule(2, N, colony, vnum);
xlabel('time[s]')
ylabel('edge index')
legend('veh 1','veh2')

%%

figure(4)
clf
% drawStl(stlAddr,4);
hold on
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
drawnow

% animate
tour = colony.queen.tour;
tick = colony.queen.tickHistory;
reservation = colony.queen.reservation;
occupancy = colony.queen.occupancy;
tourLen = colony.queen.vehTourLen;
finished = zeros(vnum,1);

% vid = VideoWriter('animation.avi','MPEG-4');
% vid.FrameRate = fps;
% vid.Quality = 100;
vid = VideoWriter('animation.avi');
open(vid);

T = max(tick,[],'all');
stepNum = ceil(T/simStep);
simT = linspace(0,T,stepNum);

for i = 1:vnum
    vv(i) = plot3(node(1,1),node(1,2),node(1,3),'o','LineWidth',3);
    hold on
end
hh(i) = line([0,0],[0,0],[0,0]);
view(0,90)
% view(-80,-10)
tt = text(0,120,num2str(0.0));

count = 0;
for t = simT
    count = count + 1;
    for v = 1:vnum
        % calc status
        routeIdx = getRouteStep(tick(v,:),t);
        if routeIdx ~= -1
            initNode = tour(v,routeIdx); termNode = tour(v,routeIdx+1);
            if A(initNode, termNode) ~= 0 && termNode ~= 1
                actualStart = occupancy{v,routeIdx+1}(1,3);
                regionDur = tick(v,routeIdx+1) - actualStart;
                weight = (t-actualStart)/regionDur;
                if t-actualStart > regionDur
                    weight = 1;
                elseif t-actualStart < 0
                    weight = 0;
                end
                newPos = (1-weight) .* node(initNode,:) + weight .* node(termNode,:);
                set(vv(v), 'XData', newPos(1), 'YData', newPos(2), 'ZData', newPos(3));
            elseif termNode == 1
                finished(v) = true;
                set(vv(v), 'XData', node(tour(v,tourLen(v)),1), 'YData', node(tour(v,tourLen(v)),2), 'ZData', node(tour(v,tourLen(v)),3),'Color',[0.4 0.4 0.4]);
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
            set(vv(v),'Color',[0.4 0.4 0.4])
        else
            set(vv(v), 'XData', node(1,1), 'YData', node(1,2), 'ZData', node(1,3));
        end
    end
    set(tt,'String',num2str(t))
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
    if finished == 1
        break;
    end
    drawnow;
    frame = getframe(gcf);
    writeVideo(vid,frame);
    pause(0.01)
end
close(vid);

%%

drawBestTour_forSoleVRPCF(colony,mapGraph,vnum,1)

figure(1)
clf
hold on
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
title('Result')
for i = 1:N
    text(node(i,1)+0.1,node(i,2),num2str(i))
end
% view(0, 90)