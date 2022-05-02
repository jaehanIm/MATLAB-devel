function superRoute = HLP_solver_ACS(map,antNo,termCond)

mapGraph = map;
N = map.n;
vnum = map.vnum;

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

bestFitness = inf;
bestTour = [];
history = zeros(maxIter,1);
tic
count = 0;
disp('Initiating HLP ACS!');

for t = 1 : maxIter
    % Create Ants 
    
    colony = [];
    colony = createColonyACSVRP( mapGraph, colony , antNo, tau, eta, alpha,  beta, q, psi, vnum);
    
    
    % Calculate the fitness values of all ants 
    for i = 1 : antNo 
        colony.ant(i).fitness = fitnessFunctionVRP_ND(colony.ant(i).tour, colony.ant(i).vehTourLen, mapGraph, vnum, map.ND);
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
    tau = updatePhromoneACSVRP(tau , colony, rho, minIndex, vnum);  
    
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
    if count >= termCond
        disp('Solution stabilized. Terminating HLP ACS!');
        break;
    end
end

superRoute = bestTour;
for i = 1:vnum
    deleteHomeIdx = find(superRoute(i,:) == 1);
    superRoute(i,deleteHomeIdx(end)) = 0;
end
superRoute(:,end) = [];

% soleACO_time = toc;
% soleACO_result = colony.queen.fitness;
% historyACO = history;
% historyACO(historyACO==0) = [];



end