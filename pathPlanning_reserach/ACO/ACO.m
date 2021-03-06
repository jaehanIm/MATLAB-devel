clear all
% close all
clc

%% Problem preparation 

% Create the graph 
[ graph ]  = createGraph();

% Draw the graph 
% figure 
 
% subplot(1,3,1)
% drawGraph( graph); 


%% ACO algorithm 
global exploreP
%% Initial parameters of ACO 
maxIter = 500;
antNo = 50;

tau0 = 10 * 1 / (  graph.n * mean( graph.edges(:)  )  );  % Initial phromone concentration

tau = tau0 * ones( graph.n , graph.n); % Phromone matirx 
eta = 1./ graph.edges;  % desirability of each edge 
% eta = ones(graph.n,graph.n);

rho = 0.5; % Evaporation rate 
alpha = 1;  % Phromone exponential parameters 
beta = 1;  % Desirability exponetial paramter
exploreP = 0.05;


%% Main loop of ACO 

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);
for t = 1 : maxIter
    % Create Ants 
    
    colony = [];
    colony = createColony( graph, colony , antNo, tau, eta, alpha,  beta);
    
    
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        colony.ant(i).fitness = fitnessFunction(colony.ant(i).tour , graph );
    end
    
    % Find the best ant (queen)
    allAntsFitness = [ colony.ant(:).fitness ];
    [ minVal , minIndex ] = min( allAntsFitness );
    if minVal < bestFitness 
        bestFitness = colony.ant(minIndex).fitness;
        bestTour = colony.ant(minIndex).tour;
    end
    
    colony.queen.tour = bestTour;
    colony.queen.fitness = bestFitness;
        
    % Update phromone matrix 
    tau = updatePhromone( tau , colony );  
    
    % Evaporation 
    tau  = ( 1 - rho ) .* tau;
    
    % Display the results 
    
    outmsg = [ 'Iteration #' , num2str(t) , ' Shortest length = ' , num2str(colony.queen.fitness)  ];
    disp(outmsg)
    history(t,1) = colony.queen.fitness;
%     subplot(1,3,1)
%     title(['Iteration #' , num2str(t) ])
%     % Visualize best tour and phromone concentration
%     subplot(1,3,2)
%     cla
%     drawBestTour( colony, graph );
%     
%     
%     subplot(1,3,3)
%     cla
%     drawPhromone( tau , graph );
%    
%    drawnow
end

drawBestTour2( colony, graph );

figure(10)
hold on
grid on
plot(history,'x--')







