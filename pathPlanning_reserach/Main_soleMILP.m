addpath('./ACO')
addpath('./MILP')
addpath('./ComDetTBv090/')
addpath('./ComDetTBv090/Algorithms/')
addpath('./ComDetTBv090/Auxiliary/')

%% Param setting
vnum = 3;
capacity = 1e10;

%% wp generator
poc_path_planner

node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];

N = size(node,1);
A = zeros(N,N); % connectivity matrix
C = zeros(N,N); % cost matrix
L = zeros(N,N); % linear distance matrix

for i = 1:N
    for j = 1:N
        L(i,j) = distance(node(i,:), node(j,:));
    end
end

%% Build Network Graph
for i = 1:N % do not connect depot!
    for j = 1:N
        if i~=j
            if L(i,j) < conThres
                A(i,j) = 1;
                C(i,j) = L(i,j);
            end
        else
            A(i,j) = 0;
        end
    end
end


G = graph(C);
degree = centrality(G,'degree');
closeness = centrality(G,'closeness');
betweenness = centrality(G,'betweenness');
pagerank = centrality(G,'pagerank');
eigenvector = centrality(G,'eigenvector');

figure(2)
p = plot(G,'Layout','force','EdgeAlpha',0.3,'MarkerSize',7);
p.NodeCData = betweenness;
colormap jet
colorbar

figure(3)
clf
grid on
hold on
plot(normalize(degree, 'range'))
plot(normalize(closeness, 'range'))
plot(normalize(betweenness, 'range'))
plot(normalize(pagerank, 'range'))
plot(normalize(eigenvector, 'range'))
legend('degree','close','betweenness','pagerank','eigen')
ylim([0 1.5])

%% Save Network
graph1.n = N;
graph1.node.x = node(:,1);
graph1.node.y = node(:,2);
graph1.node.z = node(:,3);
graph1.edges = L;

save('graph_complete.mat','graph1');

graph2 = graph(C);

%% Completefication
tic
implicitPath = [];
for i = 1:N-1
    for j = i+1:N
        if A(i,j) == 0 % if disconnected
            [path,d] = shortestpath(graph2,i,j);
            C(i,j) = d;
            C(j,i) = d;
            implicitPath{i,j} = path;
            implicitPath{j,i} = path;
        end
    end
end
toc

%% solve
map.vnum = vnum;
map.N = N;
map.capacity = capacity;
map.C = C;
map.oblig = [];

tic
[routeResult,score]=TSP_solver(map);
toc

fig2 = figure(11);
fig3 = copyobj(mapfig,fig2);
hold on
plot()