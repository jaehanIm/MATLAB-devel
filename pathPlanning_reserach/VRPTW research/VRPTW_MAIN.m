addpath('./../')

%%%%%%%%%%% ACO for VRPTW version %%%%%%%%%%%

fovFactor = 2;
mapheight = 3;
inpection_dist = 7;

distThres = 10;
vnum = 3;
antNo = 20;
stopThres = 100;

% generate map
mapGenerator_VRPTW
node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
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
for i = 1:N-1
    for j = 2:N
        if L(i,j) < distThres
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end

tic
C = zeros(N,N);
implicitRoute = cell(N,N);
for i = 1:N-1
    for j = 2:N
        if A(i,j) == 1
            C(i,j) = L(i,j);
        else
            [implRoute, implCost, ~] = shortestpath(graph(A),i,j);
            implicitRoute{i,j} = implRoute;
            implicitRoute{j,i} = fliplr(implRoute);
            A(i,j) = 1; A(j,i) = 1;
            C(i,j) = implCost; C(j,i) = implCost;
        end
    end
end
toc

% mapData generation
mapGraph.n = N;
mapGraph.edges = L;
mapGraph.node = node;

% time constraint construction
[~,~,timeMax] = NNHeuristic_VRPTW(mapGraph);


% solve problem
ACS_VRPTW;
