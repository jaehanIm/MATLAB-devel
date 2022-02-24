function [tour,score] = LLP_solver(map,maxIter,antNo)
% require N C vehNum homeIdx totSubProbNodeIdx subProbEndNodeIdx

%% Parameter Setting
global debugTemp

mapGraph.n = map.N; % Total number of subnodes
mapGraph.edges = map.C;

tau0 = 10 * 1 / (  mapGraph.n * mean( mapGraph.edges(:)  )  );  % Initial pheromone concentration
tau = tau0 * ones( mapGraph.n , mapGraph.n); % Pheromone matrix
eta = 1./ mapGraph.edges;  % desirability of each edge 

rho = 0.02; % Evaporation rate 
alpha = 1;  % Pheromone exponential parameters 
beta = 1;  % Desirability exponetial paramter
homeIdx = 1;
vehNum = 1;
debugTemp = [];

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);

for t = 1 : maxIter
%     colonyIter = [];
    mapGraph.totSubProbNodeIdx = map.totSubProbNodeIdx;
    mapGraph.subProbEndNodeIdx = map.subProbEndNodeIdx;

    %% ACO Main
    % Create tours
    colony = [];
    colony = createColonyLLP(mapGraph,colony,antNo,tau,eta,alpha,beta); %result based on local idx.
    
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        colony.ant(i).fitness = fitnessFunctionLLP(colony.ant(i).tour, colony.ant(i).vehStepLen, mapGraph);
    end

    % Find the best ant (queen)
    allAntsFitness = [ colony.ant(:).fitness ];
    [ minVal , minIndex ] = min( allAntsFitness );
    if minVal < bestFitness 
        bestFitness = colony.ant(minIndex).fitness;
        bestTour = colony.ant(minIndex).tour;
        bestTourLen = colony.ant(minIndex).vehStepLen;
    end
    
    colony.queen.tour = bestTour;
    colony.queen.fitness = bestFitness;
    colony.queen.vehStepLen = bestTourLen;
        
    % Update phromone matrix 
    tau = updatePhromoneLLP(tau , colony);  
    
    % Evaporation 
    tau  = ( 1 - rho ) .* tau;

    % Display the results
    outmsg = [ 'Iteration #' , num2str(t) , ' Shortest length = ' , num2str(colony.queen.fitness)  ];
    disp(outmsg)
    history(t,1) = colony.queen.fitness;
end

tour = colony.queen.tour;
score = colony.queen.fitness;

end
