% clear all
% close all
clc

%% Problem preparation 

% Create the graph 
graph  = createGraph();


%% ACOVRP algorithm 

%% Initial parameters of ACO 
maxIter = 300;
antNo = 100;

tau0 = 10 * 1 / (  graph.n * mean( graph.edges(:)  )  );  % Initial phromone concentration

tau = tau0 * ones( graph.n , graph.n); % Phromone matirx 
eta = 1./ graph.edges;  % desirability of each edge 
% eta = ones(graph.n,graph.n);

rho = 0.02; % Evaporation rate 
alpha = 1;  % Phromone exponential parameters 
beta = 1;  % Desirability exponetial paramter

global homeIdx vehNum vehCapacity temp  mutationRate
homeIdx = 16; %60 is the main
vehNum = 4;
vehCapacity = [100 100 100 100 100];
mutationRate = 0.0;

%% Main loop of ACO 

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);
for t = 1 : maxIter
    % Create Ants 
    
    colony = [];
    colony = createColonyVRP( graph, colony , antNo, tau, eta, alpha,  beta);
    
    
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        colony.ant(i).fitness = fitnessFunctionVRP(colony.ant(i).tour, colony.ant(i).vehTourLen, graph);
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

drawBestTourVRP( colony, graph );

figure(10)
hold on
grid on
plot(history,'x--')





% bestScore = horzcat(bestScore,colony.queen.fitness);

