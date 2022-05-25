
benchTime = solveTime;
benchScore = totalScoreHistory;
benchScoreL = totalScoreHistoryL;

mapGraph = graph1;

%% ACOVRP algorithm 

global homeIdx temp tau0
%% Initial parameters of ACO 
maxIter = 100000;

tau0 = 1 / (mapGraph.n * mean( mapGraph.edges(:)) * N);  % Initial phromone concentration

tau = tau0 * ones( mapGraph.n , mapGraph.n); % Phromone matirx 
eta = 1./ mapGraph.edges;  % desirability of each edge 

rho = 0.1; % Evaporation rate
alpha = 1;  % Phromone exponential parameters 
beta = 2;   % Desirability exponetial paramter
psi = 0.1;    % local pheromone evaporation rate
q = 0.9;    % Exploration Exploitation parameter

homeIdx = 1; %60 is the main
temp = [];

%% Main loop of ACS
oncePassed1 = 0;
oncePassed2 = 0;
oncePassed3 = 0;

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);
tic
count = 0;
% disp('Initiating ACS!');

for t = 1 : maxIter
    % Create Ants 
    
    colony = [];
    colony = createColonyACSVRP( mapGraph, colony , antNo, tau, eta, alpha,  beta, q, psi, vnum);
    
    
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        [colony.ant(i).fitness, colony.ant(i).fitnessL, colony.ant(i).fitnessPer, colony.ant(i).fitnessM] = fitnessFunctionVRP(colony.ant(i).tour, colony.ant(i).vehTourLen, mapGraph, vnum);       
    end
    
    % Find the best ant (queen)
    allAntsFitness = [ colony.ant(:).fitness ];
    [ minVal , minIndex ] = min( allAntsFitness );
    if minVal < bestFitness 
        bestFitness = colony.ant(minIndex).fitness;
        bestTour = colony.ant(minIndex).tour;
        bestTourLen = colony.ant(minIndex).vehTourLen;
        bestFitnessL = colony.ant(minIndex).fitnessL;        
        bestFitnessPer = colony.ant(minIndex).fitnessPer;
    end
    
    colony.queen.tour = bestTour;
    colony.queen.fitness = bestFitness;
    colony.queen.vehTourLen = bestTourLen;
    colony.queen.fitnessL = bestFitnessL;
    colony.queen.fitnessPer = bestFitnessPer;
        
    % Update phromone matrix 
    tau = updatePhromoneACSVRP(tau , colony, rho, minIndex, vnum);  
    
    % Display the results
    if mod(t,50) == 0
%         outmsg = [ 'Iteration #' , num2str(t) , ' Shortest length = ' , num2str(colony.queen.fitness)  ];
%         disp(outmsg)
    end
    history(t) = colony.queen.fitness;
    if t ~= 1
        if colony.queen.fitness == history(t-1)
            count = count + 1;
        else
%             disp('solution update detected!');
            count = 0;
        end
    end

    if bestFitness < benchScore && oncePassed1 == 0
        equiPerformanceTime = toc;
        oncePassed1 = 1;
    end

    if bestFitnessL < benchScoreL && oncePassed2 == 0
        equiPerformanceTimeL = toc;
        oncePassed2 = 1;
    end

    if benchTime <= toc && oncePassed3 == 0
        equiTimePerformance = bestFitness;
        oncePassed3 = 1;
    end

    if count >= stopThres
%         disp('Solution stabilized. Terminating ACS!');
        break;
    end
end
soleACO_time = toc;
soleACO_result = colony.queen.fitness;
soleACO_resultL = colony.queen.fitnessL;
soleACO_resultPerV = colony.queen.fitnessPer;

history(t+1,:) = [];
historyACO = history;
historyACO(historyACO==0) = [];

