%% Param setting
vnum = 3;
% depotPos = [60 60 60];
% depotPos = [10 -5 0];
depotPos = [-30,20,-10]; % factory


fovFactor = 1.8;
inpection_dist = 1; % Inspection distance
mapheight = 5.0;
conThres = 1.7;
% stlAddr = '/home/jaehan/Desktop/generic.stl';
% stlAddr = '/home/jaehan/Downloads/generic_edited.stl';
% stlAddr = 'C:\Users\dlawo\Desktop\Model\generic.stl';

% 최신 데이터는 totladata3
%% wp generator

node = load('stlnode.mat');
node = node.node;


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
% for i = 2:N % do not connect depot!
%     for j = 2:N
%         if i~=j
%             if L(i,j) < conThres
%                 A(i,j) = 1;
%                 C(i,j) = L(i,j);
%             end
%         else
%             A(i,j) = 0;
%         end
%     end
% end
% 
% C(1,:) = L(1,:);
% C(:,1) = L(:,1);
% 
% A_orig = A;
% C_orig = C;

%for factory
for i = 2:size(extNodes,1)
    for j = 2:size(extNodes,1)
        if i~=j
            if L(i,j) <= 5%4.3
                A(i,j) = 1;
                C(i,j) = L(i,j);
            end
        else
            A(i,j) = 0;
            C(i,j) = 0;
        end
    end
end
for i = size(extNodes,1)+1:size(extNodes,1)+size(intNodes,1)
    for j = size(extNodes,1)+1:size(extNodes,1)+size(intNodes,1)
        if i~=j
            if L(i,j) <= 3 %3
                A(i,j) = 1;
                C(i,j) = L(i,j);
            end
        else
            A(i,j) = 0;
            C(i,j) = 0;
        end
    end
end
for i = size(extNodes,1)-8:size(extNodes,1)+8
    for j = size(extNodes,1)-8:size(extNodes,1)+8
        if i~=j
            if L(i,j) <= 5
                A(i,j) = 1;
                C(i,j) = L(i,j);
            end
        else
            A(i,j) = 0;
            C(i,j) = 0;
        end
    end
end

% [A,C]=graphSparseConnection(node,A,C,L);
% A(1,:) = 0; A(:,1) = 0;
% C(1,:) = 0; C(:,1) = 0;

[A_ext,C_ext]=graphSparseConnection(node(1:size(extNodes,1),:,:),A(1:size(extNodes,1),1:size(extNodes,1)),C(1:size(extNodes,1),1:size(extNodes,1)),L(1:size(extNodes,1),1:size(extNodes,1)));
[A_int,C_int]=graphSparseConnection(node(size(extNodes,1)+1:N,:,:),A(size(extNodes,1)+1:N,size(extNodes,1)+1:N),C(size(extNodes,1)+1:N,size(extNodes,1)+1:N),L(size(extNodes,1)+1:N,size(extNodes,1)+1:N));
A(1:size(extNodes,1),1:size(extNodes,1)) = A_ext;
C(1:size(extNodes,1),1:size(extNodes,1)) = C_ext;
A(size(extNodes,1)+1:N,size(extNodes,1)+1:N) = A_int;
C(size(extNodes,1)+1:N,size(extNodes,1)+1:N) = C_int;

A(1,:) = 0; A(:,1) = 0;

A(1:size(extNodes,1)-9,size(extNodes,1)+9:N) = 0; % only for factory
A(size(extNodes,1)+9:N,1:size(extNodes,1)-9) = 0; % only for factory
C(1:size(extNodes,1)-9,size(extNodes,1)+9:N) = 0; % only for factory
C(size(extNodes,1)+9:N,1:size(extNodes,1)-9) = 0; % only for factory

A_orig = A;
C_orig = C;

G = graph(C);

totalDegree = sum(A(2:end,2:end),'all')/2
completeDegree = nchoosek(N-1,2)
degreeConnectivity = totalDegree / completeDegree

% G = graph(C);
% degree = centrality(G,'degree');
% closeness = centrality(G,'closeness');
% betweenness = centrality(G,'betweenness');
% pagerank = centrality(G,'pagerank');
% eigenvector = centrality(G,'eigenvector');
% avgdegree = mean(degree);
% 
% figure(2)
% p = plot(G,'Layout','force','EdgeAlpha',0.3,'MarkerSize',7);
% p.NodeCData = betweenness;
% colormap jet
% colorbar
% 
% figure(3)
% clf
% grid on
% hold on
% plot(normalize(degree, 'range'))
% plot(normalize(closeness, 'range'))
% plot(normalize(betweenness, 'range'))
% plot(normalize(pagerank, 'range'))
% plot(normalize(eigenvector, 'range'))
% legend('degree','close','betweenness','pagerank','eigen')
% ylim([0 1.5])


%% Completefication
tic
disp("complefication start")
graph2 = graph(C);
implicitPath = [];
for i = 2:N-1
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
for i = 1:N
    C(1,i) = L(1,i);
    C(i,1) = L(1,i);
end
A(1,:) = 1; A(:,1) = 1;
disp("complefication complete")
completeTime_soleACO = toc;

%% Save Network
graph1.n = N;
graph1.node.x = node(:,1);
graph1.node.y = node(:,2);
graph1.node.z = node(:,3);
graph1.edges = C;

save('graph_complete.mat','graph1');

%% solve
antNo = 20;
stopThres = 400;
% ACOVRP_forSoleACO

ACSVRP_forSoleACS

% IPt = [0.02,0.69,5.71,21.6,197.24,1825.75];
% ACOt = [3.78,5.67,7.83,7.4,10.99,28.91,58.2,58.13,150.26,201.57];
% IPs = [91.08,131.7,180.83,127.44,170.93,157.73];
% ACOs = [91.08,131.7,180.89,127.44,173.81,163.50,206.07,304.03,345.25,458.35];
% OptGap = (ACOs(1:6)-IPs)./IPs * 100;
% IPlab = [6 9 12 12 15 20];
% ACOlab = [6 9 12 12 15 20 30 63 88 130];
% figure(1)
% clf
% yyaxis left
% loglog(IPlab,IPt,'bo-','LineWidth',2.5)
% hold on
% grid on
% loglog(ACOlab,ACOt,'ro-','LineWidth',2.5)
% xlabel('Node Number (log scale)')
% ylabel('Computation Time [s]')
% title('Node Number to Computation Time & Optimality Gap')
% yyaxis right
% loglog(IPlab,OptGap,'kx-','LineWidth',2.5,'MarkerSize',10)
% ylabel('optimality gap [%]')
% ylim([0 10])
% left_color = [0 0 0];
% right_color = [0 0 0];
% set(figure(1),'defaultAxesColorOrder',[left_color; right_color]);
% legend('Exact','ACS','Optimality Gap')

disp('Time')
completeTime_soleACO
soleACO_time
soleACO_result

drawBestTour_forSoleACO( colony, mapGraph, vnum);
drawStl(stlAddr,1)

