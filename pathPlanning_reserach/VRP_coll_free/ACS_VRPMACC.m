
%% ACOVRP algorithm 
global homeIdx tau0
stopThres = 200;
%% Initial parameters of ACO 
maxIter = 10000;

tau0 = 1 / (mapGraph.n * mean( mapGraph.edges(:)) * N);  % Initial phromone concentration
tau = tau0 * ones( mapGraph.n , mapGraph.n); % Phromone matirx
eta = 1./ mapGraph.edges;  % desirability of each edge

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
    
    colony = [];
    colony = createColonyVRPMACC( mapGraph, colony , antNo, tau, eta, param, vnum);
    
    
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        [colony.ant(i).fitness, colony.ant(i).fitnessL, colony.ant(i).fitnessPer, colony.ant(i).fitnessM] = ...
            fitnessFunctionVRPMACC(colony.ant(i).tour, colony.ant(i).vehTourLen, mapGraph, vnum);       
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
    tau = updatePhromoneVRPMACC(tau , colony, param.rho, minIndex, vnum);  
    
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
%             disp('solution update detected!');
            count = 0;
        end
    end

    if count >= stopThres
        disp('Solution stabilized. Terminating ACS!');
        break;
    end
end
soleACO_time = toc;