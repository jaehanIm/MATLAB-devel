%% Param setting
vnum = 3;
% depotPos = [60 60 60];
depotPos = [10 -5 0];

fovFactor = 1.8;
inpection_dist = 1; % Inspection distance
mapheight = 5.0;
conThres = 1.7;
stlAddr = '/home/jaehan/Desktop/generic.stl';
% stlAddr = '/home/jaehan/Downloads/generic_edited.stl';
stlAddr = 'C:\Users\dlawo\Desktop\Model\generic.stl';


%% wp generator

node = loadStl(stlAddr,inpection_dist);

% node segmentation
body = node;

tail = body(body(:,1) > 14.7421 & body(:,2) < 0.955 & body(:,2) > -0.955,:,:);
body(body(:,1) > 14.7421 & body(:,2) < 0.955 & body(:,2) > -0.955,:,:) = [];
factor = 0.1; localN = size(tail,1);
tail = tail(rand(localN,1) < factor,:,:);

lEngine = body(body(:,2) < -1.9 & body(:,2) > -5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) < 0.22,:,:);
body(body(:,2) < -1.9 & body(:,2) > -5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) < 0.22,:,:) = [];
factor = 0.1; localN = size(lEngine,1);
lEngine = lEngine(rand(localN,1) < factor,:,:);

rEngine = body(body(:,2) > 1.9 & body(:,2) < 5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) < 0.22,:,:);
body(body(:,2) > 1.9 & body(:,2) < 5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) < 0.22,:,:) = [];
factor = 0.1; localN = size(rEngine,1);
rEngine = rEngine(rand(localN,1) < factor,:,:);

emphennage = body(body(:,1)>12 & body(:,2)<2 & body(:,2)>-2,:,:);
body(body(:,1)>12 & body(:,2)<2 & body(:,2)>-2,:,:) = [];
factor = 0.2; localN = size(emphennage,1);
emphennage = emphennage(rand(localN,1) < factor,:,:);

lEngineNeg = body(body(:,2) < -1.9 & body(:,2) > -5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) > 0.22 & body(:,3) < 0.5,:,:);
body(body(:,2) < -1.9 & body(:,2) > -5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) > 0.22 & body(:,3) < 0.3,:,:) = [];

rEngineNeg = body(body(:,2) > 1.9 & body(:,2) < 5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) > 0.22 & body(:,3) < 0.5,:,:);
body(body(:,2) > 1.9 & body(:,2) < 5.13 & body(:,1) > 1.13 & body(:,1) < 4.08 & body(:,3) > 0.22 & body(:,3) < 0.3,:,:) = [];

body(body(:,2) > -3.59 & body(:,2) < -3.51 & body(:,1) > 1.5 & body(:,1) < 3.95 & body(:,3) > 0.22,:,:) = [];
body(body(:,2) > 3.51 & body(:,2) < 3.59 & body(:,1) > 1.5 & body(:,1) < 3.95 & body(:,3) > 0.22,:,:) = [];

nose = body(body(:,1) < -2.9,:,:);
body(body(:,1) < -2.9,:,:) = [];
factor = 0.5; localN = size(nose,1);
nose = nose(rand(localN,1) < factor,:,:);

node = [body;tail;lEngine;rEngine;nose;emphennage];


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

C(1,:) = L(1,:);
C(:,1) = L(:,1);

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

