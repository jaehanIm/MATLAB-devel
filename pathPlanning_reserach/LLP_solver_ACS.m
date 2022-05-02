function [tour,score, clusterCost, bridgeCost, residueCost] = LLP_solver_ACS(map,antNo,termCond)
% require N C vehNum homeIdx totSubProbNodeIdx subProbEndNodeIdx

%% Parameter Setting
global debugTemp tau0

mapGraph.n = map.N; % Total number of subnodes
mapGraph.edges = map.C;

tau0 = 1 * 1 / (  mapGraph.n * mean( mapGraph.edges(:)  )  );  % Initial pheromone concentration
tau = tau0 * ones( mapGraph.n , mapGraph.n); % Pheromone matrix
eta = 1./ mapGraph.edges;  % desirability of each edge 

rho = 0.1; % Evaporation rate
alpha = 1;  % Phromone exponential parameters 
beta = 2;   % Desirability exponetial paramter
psi = 0.1;    % local pheromone evaporation rate
q = 0.9;    % Exploration Exploitation parameter
homeIdx = 1;
vnum = 1;
debugTemp = [];

bestFitness = inf;
bestTour = [];
count = 1;
maxIter = 1e4;
history = zeros(maxIter,1);

disp('Initiating ACS!');
for t = 1 : maxIter
%     colonyIter = [];
    mapGraph.totSubProbNodeIdx = map.totSubProbNodeIdx;
    mapGraph.subProbEndNodeIdx = map.subProbEndNodeIdx;

    %% ACO Main
    % Create tours
    colony = [];
    colony = createColonyACSLLP(mapGraph,colony,antNo,tau,eta,alpha,beta,q,psi); %result based on local idx.
    
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
        bestBridgeCost = colony.ant(minIndex).bridgeCost;
        bestClusterCost = colony.ant(minIndex).clusterCost;
    end
    
    colony.queen.tour = bestTour;
    colony.queen.fitness = bestFitness;
    colony.queen.vehStepLen = bestTourLen;
    colony.queen.bridgeCost = bestBridgeCost;
    colony.queen.clusterCost = bestClusterCost;
        
    % Update phromone matrix 
    tau = updatePhromoneACSLLP(tau , colony, rho, minIndex);  
    
    % Evaporation 
    tau  = ( 1 - rho ) .* tau;

    % Display the results
    if mod(t,50) == 0
        outmsg = [ 'Iteration #' , num2str(t) , ' Shortest length = ' , num2str(colony.queen.fitness)  ];
        disp(outmsg)
    end

    history(t,1) = colony.queen.fitness;

    if t ~= 1
        if colony.queen.fitness == history(t-1)
            count = count + 1;
        else
            disp('solution update detected!');
            count = 0;
        end
    end
    if count >= termCond
        disp('Solution stabilized. Terminating ACS!');
        disp('=====================================');
        break;
    end

end

%% Scoring

history(count+1:end,:) = [];
tour = colony.queen.tour;

initCost = map.C(tour(1),tour(2));
returnCost = map.C(tour(end-1),tour(end));
colony.queen.bridgeCost(1) = [];

score = [colony.queen.fitness];
bridgeCost = colony.queen.bridgeCost;
clusterCost = colony.queen.clusterCost;
residueCost = [initCost,returnCost];

end
