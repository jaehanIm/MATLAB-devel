%% Param setting
vnum = 3;
% depotPos = [60 60 60];
depotPos = [10 -5 0];

fovFactor = 1.8;
inpection_dist = 1; % Inspection distance
mapheight = 5.0;
conThres = 3;
stlAddr = '/home/jaehan/Desktop/generic.stl';
% stlAddr = '/home/jaehan/Downloads/generic_edited.stl';

%% wp generator

node = loadStl(stlAddr,inpection_dist);
N = size(node,1);

% node_temp = zeros(floor(N/2),3);
% for i = 1:floor(size(node_temp,1))
%     node_temp(i,:) = node(i*2-1,:);
% end
% node = node_temp;
% node(find(isnan(node(:,1))),:) = [];

% node reduction
N = size(node,1);
factor = 2;
node_temp = zeros(floor(N/factor),3);
for i = 1:floor(size(node_temp,1))
    node_temp(i,:) = node(i*factor-1,:);
end
node = node_temp;
node(find(isnan(node(:,1))),:) = [];

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
for i = 2:N % do not connect depot!
    for j = 2:N
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

A_orig = A;
C_orig = C;


G = graph(C);
degree = centrality(G,'degree');
closeness = centrality(G,'closeness');
betweenness = centrality(G,'betweenness');
pagerank = centrality(G,'pagerank');
eigenvector = centrality(G,'eigenvector');
avgdegree = mean(degree);

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


%% Completefication
tic
disp("complefication start")
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
disp("complefication complete")
completeTime_soleACO = toc

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

