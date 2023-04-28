
addpath('./../')
addpath('..\ACO\')
addpath('./../Model/')

%%%%%%%%%%% ACO for VRPCF version %%%%%%%%%%%

% parameter setting
fovFactor = 2.7;
mapheight = 3;
inpection_dist = 7;

antNo = 20;
stopThres = 20;
capacity = 395;
servTime = 1;
homePos = [0,0,0];


%% random node generator
N = 40;
node = rand(N,3);
node(:,1:2) = node(:,1:2) * 100;
node(:,3) = node(:,3) * 10;
node = vertcat(homePos,node);
simStep = 1;
distThres = 60;
vnum =6;

% mapGenerator_VRPCF
% node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
% node = vertcat(homePos,node);
% simStep = 1;
% distThres = 20;
% vnum = 6;

% node = [1,0,0;2,0,0;3,0,0;4,0,0;2,1,0];
% node = vertcat(homePos,node);
% simStep = 0.1;
% distThres = 20;
% vnum = 1;

N = size(node,1);

%% Graph construction
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
G = graph(C);
tic
implicitRoute = cell(N,N);
for i = 1:N-1
    for j = i+1:N
        if A(i,j) == 0
            [implRoute, implCost, ~] = shortestpath(G,i,j);
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

%% Solve
ACS_VRPDM;
colony.queen.fitnessL
colony.queen.fitnessDev

%% Draw result
solution = colony.queen.tour;
detailRoute = cell(vnum,1);
for v = 1:vnum
    detailRoute{v} = routeDecoder(nonzeros(solution(v,:)), implicitRoute,N);
end
arrivTime = zeros(vnum,1);
for i = 1:vnum
    arrivTime(i) = getTimeAtArrival(solution(i,:), C);
end
maxT = max(arrivTime);
simT = 0:simStep:maxT;

% draw map
figure(4)
clf
hold on
for i = 1:size(node,1)
    plot3(node(:,1),node(:,2),node(:,3),'.','MarkerSize',10,'Color',"#A2142F");
end
view(0, 90)
lineList = find(A(:));
for i = 1:N-1
    for j= i+1:N
        if A(i,j) ~= 0
            initIdx = i;
            termIdx = j;
            line([node(initIdx,1), node(termIdx,1)],[node(initIdx,2), node(termIdx,2)],[node(initIdx,3), node(termIdx,3)]-10,'Color',[0.7 0.7 0.7]);
        end
    end
end
hold on
plot3(node(1,1),node(1,2),node(1,3),'x','MarkerSize',5,'LineWidth',4)
grid on
axis equal
title('Tour animation')
view(0, 90)
drawnow

% open video
vid = VideoWriter('forfun11.avi');
open(vid);

% draw initial
for v = 1:vnum
    vv(v) = plot3(node(1,1),node(1,2),node(1,3),'o','LineWidth',3);    
end
tt = text(homePos(1)+2,homePos(2),num2str(0.0));

% update state
for t = simT
    for v = 1:vnum
        if t < arrivTime(v)
            pos = getVehPos(detailRoute{v}, t, C, node, implicitRoute, N);
        else
            pos = [homePos(1),homePos(2),homePos(3)];
        end
        set(vv(v),'XData',pos(1),'YData',pos(2),'ZData',pos(3));
        set(tt,'String',num2str(t))
    end
    frame = getframe(gcf);
    writeVideo(vid,frame);
    drawnow;
end
close(vid)