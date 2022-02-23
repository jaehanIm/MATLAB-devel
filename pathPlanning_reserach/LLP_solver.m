function LLP_solver(mapGraph,maxIter,antNo)
% require N C vehNum homeIdx endIdx

%% ACOVRP algorithm

mapGraph.n = mapGraph.N;
mapGraph.edges = mapGraph.C;
mapGraph.endIdx;

%% Initial parameters of ACO 

tau0 = 10 * 1 / (  mapGraph.n * mean( mapGraph.edges(:)  )  );  % Initial pheromone concentration

tau = tau0 * ones( mapGraph.n , mapGraph.n); % Phromone matirx 
eta = 1./ mapGraph.edges;  % desirability of each edge 

rho = 0.02; % Evaporation rate 
alpha = 1;  % Phromone exponential parameters 
beta = 1;  % Desirability exponetial paramter

global debugTemp
homeIdx = mapGraph.homeIdx;
vehNum = mapGraph.vehNum;
debugTemp = [];

%% Main loop of ACO 

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);
for t = 1 : maxIter
    % Create Ants 
    colony = [];
    colony = createColonyLLP(mapGraph,colony,antNo,tau,eta,alpha,beta);
    
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        colony.ant(i).fitness = fitnessFunctionLLP(colony.ant(i).tour, colony.ant(i).vehTourLen, mapGraph);
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
    tau = updatePhromoneLLP(tau , colony, vehNum);  
    
    % Evaporation 
    tau  = ( 1 - rho ) .* tau;
    
    % Display the results
    outmsg = [ 'Iteration #' , num2str(t) , ' Shortest length = ' , num2str(colony.queen.fitness)  ];
    disp(outmsg)
    history(t,1) = colony.queen.fitness;
end

end
