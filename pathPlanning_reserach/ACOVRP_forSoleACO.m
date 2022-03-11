
graph1 = load('graph_complete.mat');
mapGraph = graph1.graph1;

% mapGraph.graph1 = G;
% mapGraph.n = N;
% mapGraph.node.x = node(:,1);
% mapGraph.node.x = node(:,2);
% mapGraph.node.x = node(:,3);
% mapGraph.edges = L;

%% ACOVRP algorithm 

%% Initial parameters of ACO 
maxIter = 10000;
antNo = 200;

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
count = 0;
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
    if mod(t,10) == 0
        outmsg = [ 'Iteration #' , num2str(t) , ' Shortest length = ' , num2str(colony.queen.fitness)  ];
        disp(outmsg)
    end
    history(t) = colony.queen.fitness;
    if t ~= 1
        if colony.queen.fitness == history(t-1)
            count = count + 1;
        else
            disp('solution update detected!');
            count = 0;
        end
    end
    if count >= 200
        disp('Solution stabilized. Terminating ACS!');
        break;
    end
end
toc
drawBestTour_forSoleACO( colony, mapGraph );
