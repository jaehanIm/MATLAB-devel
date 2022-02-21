addpath('./ACO')
addpath('./MILP')
addpath('./ComDetTBv090/')
addpath('./ComDetTBv090/Algorithms/')
addpath('./ComDetTBv090/Auxiliary/')

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
conThres = 10;
for i = 1:N
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

% figure(3)
% clf
% subplot(5,1,1)
% grid on
% plot(degree)
% ylabel('degree')
% subplot(5,1,2)
% grid on
% plot(closeness)
% ylabel('closeness')
% subplot(5,1,3)
% grid on
% plot(betweenness)
% ylabel('betweenness')
% subplot(5,1,4)
% grid on
% plot(pagerank)
% ylabel('pagerank')
% subplot(5,1,5)
% grid on
% plot(eigenvector)
% ylabel('eigenvector')

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

figure(4)
clf
grid on
hold on
p = plot3(airPosX(:),airPosY(:),airPosZ(:),'kx');
mesh(voxelPosX,voxelPosY,voxelFilterData);

for i = 1:size(G.Edges,1)
    startIdx = G.Edges.EndNodes(i,1);
    EndIdx = G.Edges.EndNodes(i,2);
    startPos = node(startIdx,:);
    EndPos = node(EndIdx,:);
    line([startPos(1) EndPos(1)],[startPos(2) EndPos(2)],[startPos(3) EndPos(3)]);
end
axis equal

%% Save Network
graph1.n = N;
graph1.node.x = node(:,1);
graph1.node.y = node(:,2);
graph1.node.z = node(:,3);
graph1.edges = L;

save('graph_complete.mat','graph1');

%% Graph Clustering

cluIdx = GCModulMax1(A);
cluNum = max(cluIdx);
nodeInCluIdx = [];
for i = 1:cluNum
    nodeInCluIdx{i} = find(cluIdx == i);
end
intraCluIdxSet = [];
intraCluIdx = [];
for i = 1:cluNum-1
    for j = i+1:cluNum
        [I,J]=find(A(nodeInCluIdx{i},nodeInCluIdx{j}));
        intraCluIdxSet{i,j} = [I,J];
        intraCluIdxSet{j,i} = [J,I];
    end
end

figure(5)
p = plot(G,'Layout','force','EdgeAlpha',0.3,'MarkerSize',7);
p.NodeCData = betweenness;
colormap jet
colorbar

figure(6)
clf
p = plot(graph(A),'Layout','force','EdgeAlpha',0.3,'MarkerSize',7);
p.NodeCData = cluIdx;
colormap jet 
colorbar

%% complete-fy

%Intra cluster complefication
