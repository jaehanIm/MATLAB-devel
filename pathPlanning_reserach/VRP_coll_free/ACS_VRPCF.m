
%% ACOVRP algorithm 

global homeIdx tau0
%% Initial parameters of ACO 
maxIter = 10000;

tau0 = 1 / (mapGraph.n * mean( mapGraph.edges(:)) * N);  % Initial phromone concentration
tau = tau0 * ones( mapGraph.n , mapGraph.n); % Phromone matirx
eta = 1./ mapGraph.edges;  % desirability of each edge
reservation = sparse(mapGraph.n, mapGraph.n);

rho = 0.1; % Evaporation rate
alpha = 1;  % Phromone exponential parameters 
beta = 1;   % Desirability exponetial paramter
gamma = 3; % time importance parameter
lambda = 1; % time scale factor
psi = 0.1;    % local pheromone evaporation rate
q = 0.9;    % Exploration Exploitation parameter
penalFactor = 2; % time violation penalty factor

homeIdx = 1; %60 is the main

%% Main loop of ACS

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);
tic
count = 0;
disp('Initiating ACS!');
colony = [];
colony.reservation = reservation;

for t = 1 : maxIter
    % Create Ants 
    colony = createColonyVRPCF( mapGraph, colony , antNo, tau, eta, alpha,  beta, gamma,...
        lambda, q, psi, vnum, timeWindow, capacity, servTime);
        
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        [colony.ant(i).fitness, colony.ant(i).fitnessL, colony.ant(i).fitnessPer, colony.ant(i).fitnessM, colony.ant(i).violation] = ...
            fitnessFunctionVRPCF(colony.ant(i).tour, colony.ant(i).vehTourLen, mapGraph, vnum, timeWindow, colony.ant(i).tickHistory, servTime, penalFactor);
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
        bestTick = colony.ant(minIndex).tick;
        bestTickHistory = colony.ant(minIndex).tickHistory;
        bestViolation = colony.ant(minIndex).violation;
    end
    
    colony.queen.tour = bestTour;
    colony.queen.fitness = bestFitness;
    colony.queen.vehTourLen = bestTourLen;
    colony.queen.fitnessL = bestFitnessL;
    colony.queen.fitnessPer = bestFitnessPer;
    colony.queen.tick = bestTick;
    colony.queen.tickHistory = bestTickHistory;
    colony.queen.violation = bestViolation;
        
    % Update phromone matrix 
    tau = updatePhromoneVRPCF(tau , colony, rho, minIndex, vnum);  
    
    % Display the results
    if mod(t,50) == 0
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

    if count >= stopThres
        disp('Solution stabilized. Terminating ACS!');
        break;
    end

end
VRPCF_solveTime = toc;
VRPCF_result = colony.queen.fitness;
VRPCF_resultL = colony.queen.fitnessL;
VRPCF_resultPerV = colony.queen.fitnessPer;

history(t+1,:) = [];
historyACO = history;
historyACO(historyACO==0) = [];

