
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

% homePos = [0,0,0];

global wowCount;
wowCount = 0;

%% random node generator
N = 80;
node = rand(N,3);
node(:,1:2) = node(:,1:2) * 100;
node(:,3) = node(:,3) * 10;
node = vertcat(homePos,node);
simStep = 1;
distThres = 20;
vnum =8;

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
ACS_VRPCF;

%% Draw result