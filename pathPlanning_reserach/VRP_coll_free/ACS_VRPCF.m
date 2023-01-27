
%% ACOVRP algorithm 
global homeIdx tau0
%% Initial parameters of ACO 
maxIter = 10000;

tau0 = 1 / (mapGraph.n * mean( mapGraph.edges(:)) * N);  % Initial phromone concentration
tau = tau0 * ones( mapGraph.n , mapGraph.n); % Phromone matirx
eta = 1./ mapGraph.edges;  % desirability of each edge

% param.rho = 0.1; % Evaporation rate
% param.alpha = 1;  % Phromone exponential parameters 
% param.beta = 1;   % Desirability exponetial paramter
% param.gamma = 3; % time importance parameter
% param.lambda = 1; % time scale factor
% param.psi = 0.1;    % local pheromone evaporation rate
% param.q = 0.9;    % Exploration Exploitation parameter
% param.penalFactor = 2; % time violation penalty factor

param.rho = 0.1; % Evaporation rate
param.alpha = 1;  % Phromone exponential parameters 
param.beta = 1;   % Desirability exponetial paramter
param.psi = 0.1;    % local pheromone evaporation rate
param.q = 0.9;    % Exploration Exploitation parameter

homeIdx = 1; %60 is the main

%% Main loop of ACS

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);
tic
count = 0;
disp('Initiating ACS!');
colony = [];
colony.reservation = [];

for t = 1 : maxIter
    % Create Ants 
    [colony, incompleteFlag] = createColonyVRPCF( mapGraph, antNo, tau, eta, param, vnum, capacity, servTime, implicitRoute);
        
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        [colony.ant(i).fitness, colony.ant(i).fitnessL, colony.ant(i).fitnessPer, colony.ant(i).fitnessM, colony.ant(i).violation] = ...
            fitnessFunctionVRPCF(colony.ant(i).tour, colony.ant(i).vehTourLen, mapGraph, vnum, colony.ant(i).tickHistory, servTime, param, incompleteFlag, N);
    end
    
    % Find the best ant (queen)
    allAntsFitness = [ colony.ant(:).fitness ];
    [ minVal , minIndex ] = min( allAntsFitness );
    if minVal < bestFitness 
        bestFitness = colony.ant(minIndex).fitness;
        bestTour = colony.ant(minIndex).tour;
        bestTourLen = colony.ant(minIndex).vehTourLen;
        bestFitnessL = colony.ant(minIndex).fitnessL;        
        bestFitnessM = colony.ant(minIndex).fitnessM;
        bestFitnessPer = colony.ant(minIndex).fitnessPer;
        bestTick = colony.ant(minIndex).tick;
        bestTickHistory = colony.ant(minIndex).tickHistory;
        bestViolation = colony.ant(minIndex).violation;
        bestReservation = colony.ant(minIndex).reservation;
        bestOccupancy = colony.ant(minIndex).occupancy;
    end
    
    colony.queen.tour = bestTour;
    colony.queen.fitness = bestFitness;
    colony.queen.vehTourLen = bestTourLen;
    colony.queen.fitnessL = bestFitnessL;
    colony.queen.fitnessM = bestFitnessM;
    colony.queen.fitnessPer = bestFitnessPer;
    colony.queen.tick = bestTick;
    colony.queen.tickHistory = bestTickHistory;
    colony.queen.violation = bestViolation;
    colony.queen.reservation = bestReservation;
    colony.queen.occupancy = bestOccupancy;
        
    % Update phromone matrix 
    tau = updatePhromoneVRPCF(tau , colony, param.rho, minIndex, vnum);
    
    % Display the results
    if mod(t,1) == 0
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

    coveredNode = unique(bestTour);
    coveredNode(1) = [];
    missingNode = ones(N,1);
    missingNode(coveredNode) = 0;
    missingNode = find(missingNode);
    if ~isempty(missingNode)
        missingNode
    end

    if count >= stopThres
        disp('Solution stabilized. Terminating ACS!');
        if ~isempty(missingNode)
            disp('Solution incomplete');
            missingNode
        else
            disp('Solution complete')
        end
        break;
    end

end
VRPCF_solveTime = toc;
VRPCF_result = colony.queen.fitness;
VRPCF_resultL = colony.queen.fitnessL;
VRPCF_resultPerV = colony.queen.fitnessPer;

