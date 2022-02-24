% function LLP_solver(map,maxIter,antNo)
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
    endNodes = map.subProbEndNodeIdx;

    %% ACO Main
    % Create tours
    colony = [];
    colony = createColonyLLP(mapGraph,colony,antNo,tau,eta,alpha,beta); %result based on local idx.

    % translate colony result
    for i = 1:antNo
        map.totSubProbNodeIdx{c}(colony.ant(i).tour);
    end

    % add colonyIter
    if c ~= 1
        for i = 1:antNo
            colonyIter.ant(i).tour = vertcat(colonyIter.ant(i).tour,map.totSubProbNodeIdx{c}(colony.ant(i).tour));
            colonyIter.ant(i).vehTourLen = colonyIter.ant(i).vehTourLen + colony.ant(i).vehTourLen;
        end
    else
        colonyIter = colony;
        for i = 1:antNo
            colonyIter.ant(i).tour = map.totSubProbNodeIdx{c}(colony.ant(i).tour);
        end
    end
end

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



% end
