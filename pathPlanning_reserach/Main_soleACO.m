addpath('./ACO')
addpath('./MILP')
addpath('./ComDetTBv090/')
addpath('./ComDetTBv090/Algorithms/')
addpath('./ComDetTBv090/Auxiliary/')

%% Param setting
% vnum = 4;
% depotPos = [60 60 60];
% antNo = 20;
% stopThres = 400;

%% wp generator
% poc_path_planner
node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
node = vertcat(depotPos,node); % add depot


% node = readtable('/home/jaehan/Desktop/TSP_test.csv');
% node = node.Variables;
% node(:,1:2) = node(:,2:3);
% node(:,3) = 0;

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

totalDegree = sum(A(2:end,2:end),'all')/2
completeDegree = nchoosek(N-1,2)
degreeConnectivity = totalDegree / completeDegree

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
completeTime_soleACO = toc

%% Save Network
graph1.n = N;
graph1.node.x = node(:,1);
graph1.node.y = node(:,2);
graph1.node.z = node(:,3);
graph1.edges = C;

save('graph_complete.mat','graph1');

%% solve
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
drawBestTour_forSoleACO( colony, mapGraph, vnum);

disp('Time')
completeTime_soleACO
soleACO_time
soleACO_result
soleACO_resultL
soleACO_resultPerV
