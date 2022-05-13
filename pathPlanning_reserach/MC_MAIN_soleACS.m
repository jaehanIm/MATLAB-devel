
%% wp generator
node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
node = vertcat(depotPos,node); % add depot

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

[A,C]=graphSparseConnection(node,A,C,L);
A(1,:) = 0; A(:,1) = 0;

G = graph(C);
A_orig = A;
C_orig = C;

%% Completefication
tic
graph2 = graph(C);
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
completeTime_soleACO = toc;

%% Save Network
graph1.n = N;
graph1.node.x = node(:,1);
graph1.node.y = node(:,2);
graph1.node.z = node(:,3);
graph1.edges = C;

% save('graph_complete.mat','graph1');

%% solve
ACSVRP_forSoleACS