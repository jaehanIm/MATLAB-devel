
% graph1 = load('graph_complete.mat');
% mapGraph = graph1.graph1;

mapGraph.graph1 = G;
mapGraph.n = N;
mapGraph.node.x = node(:,1);
mapGraph.node.x = node(:,2);
mapGraph.node.x = node(:,3);
mapGraph.edges = L;

%% ACOVRP algorithm 

%% Initial parameters of ACO 
maxIter = 150;
antNo = 50;

tau0 = 10 * 1 / (  mapGraph.n * mean( mapGraph.edges(:)  )  );  % Initial phromone concentration

tau = tau0 * ones( mapGraph.n , mapGraph.n); % Phromone matirx 
eta = 1./ mapGraph.edges;  % desirability of each edge 

rho = 0.02; % Evaporation rate 
alpha = 1;  % Phromone exponential parameters 
beta = 1;  % Desirability exponetial paramter

global homeIdx vehNum temp  mutationRate
homeIdx = 1; %60 is the main
vehNum = 3;
mutationRate = 0.0;
temp = [];

%% Main loop of ACO 

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);
tic
for t = 1 : maxIter
    % Create Ants 
    
    colony = [];
    colony = createColonyVRP( mapGraph, colony , antNo, tau, eta, alpha,  beta);
    
    
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        colony.ant(i).fitness = fitnessFunctionVRP(colony.ant(i).tour, colony.ant(i).vehTourLen, mapGraph);
    end
    
    % Find the best ant (queen)
    allAntsFitness = [ colony.ant(:).fitness ];
    [ minVal , minIndex ] = min( allAntsFitness );
    if minVal < bestFitness 
        bestFitness = colony.ant(minIndex).fitness;
        bestTour = colony.ant(minIndex).tour;
        bestTourLen = colony.ant(minIndex).vehTourLen;
    end
    
    colony.queen.tour = bestTour;
    colony.queen.fitness = bestFitness;
    colony.queen.vehTourLen = bestTourLen;
        
    % Update phromone matrix 
    tau = updatePhromoneVRP(tau , colony);  
    
    % Evaporation 
    tau  = ( 1 - rho ) .* tau;
    
    % Display the results
    outmsg = [ 'Iteration #' , num2str(t) , ' Shortest length = ' , num2str(colony.queen.fitness)  ];
    disp(outmsg)
    history(t,1) = colony.queen.fitness;
end
toc
% drawBestTourVRP( colony, mapGraph );

figure(444)
clf
grid on
hold on
for i = 1:size(G.Edges,1)
    startIdx = G.Edges.EndNodes(i,1);
    EndIdx = G.Edges.EndNodes(i,2);
    startPos = node(startIdx,:);
    EndPos = node(EndIdx,:);
    line([startPos(1) EndPos(1)],[startPos(2) EndPos(2)],[startPos(3) EndPos(3)]);
end
for i = 1:cluNum
    temp = node(nodeInCluIdx{i},:);
    plot3(temp(:,1),temp(:,2),temp(:,3),'x','LineWidth',5,'MarkerSize',5);
end
mesh(voxelPosX,voxelPosY,voxelFilterData);
plot3(superNet.pos(:,1),superNet.pos(:,2),superNet.pos(:,3),'yx','MarkerSize',10,'LineWidth',5)
for i = 1:cluNum
    text(superNet.pos(i,1),superNet.pos(i,2),superNet.pos(i,3),num2str(i));
end
colony.queen.tour(colony.queen.tour == 0) = 1;
for i = 1:vnum
    plot3(node(colony.queen.tour(i,:),1),node(colony.queen.tour(i,:),2),node(colony.queen.tour(i,:),3)+2,'r:','LineWidth',1.3)
end
axis equal

