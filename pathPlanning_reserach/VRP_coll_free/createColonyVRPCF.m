function [ colony ] = createColonyVRPCF( graph, antNo, tau, eta, param, vehNum, capacity, servTime, implicitRoute)

global homeIdx debugTemp tau0

alpha = param.alpha;
beta = param.beta;
gamma = param.gamma;
q = param.q;
lambda = param.lambda;
psi = param.psi;

A = graph.A;
N = graph.n;
C = graph.edges;
debugTemp.N = N;
colony = [];

for i = 1 : antNo
    
    unvisitedNum = N;
    
    initial_node = homeIdx;
    colony.ant(i).tour(:,1) = ones(vehNum,1) * initial_node;
    colony.ant(i).tick = zeros(vehNum,1);
    colony.ant(i).tickHistory = zeros(vehNum,1);
    unvisitedNum = unvisitedNum - 1;
    vehTourLen = ones(vehNum,1);
    reservation{N,N} = [];
    for m = 1:N
        for n = 1:N
            reservation{m,n}.num = 0;
            reservation{m,n}.info = [];
        end
    end

    while unvisitedNum ~= 0
%         for j = 1:vehNum
%             j = leastTourAgent(colony.ant(i).tour,C,vehNum,vehTourLen);
            j = getActiveVeh(colony.ant(i).tick, capacity);

            if unvisitedNum == 0
                break;
            end            
            
            % i = antNo, j = vehNum
            currentNode = colony.ant(i).tour(j,vehTourLen(j));
            P_allNodes = tau(currentNode,:).^alpha.*eta(currentNode,:).^beta;
            P_allNodes(nonzeros(colony.ant(i).tour(:))) = 0;
            P = P_allNodes./sum(P_allNodes);
            
            % for debugging
            debugTemp.P_allNodes = P_allNodes;
            debugTemp.P = P;
            debugTemp.currentNode = currentNode;
            debugTemp.i = i;
            debugTemp.j = j;
            debugTemp.colony = colony;
            debugTemp.vehTourLen = vehTourLen;
            debugTemp.unvisitedNum = unvisitedNum;
            if unvisitedNum == N-1
                debugTemp.initColony = colony;
            end
            
            if rand(1) > q
                nextNode = rouletteWheel(P);
            else
                [~,nextNode] = max(P_allNodes);
                nextNode = nextNode(1); % just in case multiple nextNodes are generated
            end

            vehTourLen(j) = vehTourLen(j) +1;
            unvisitedNum = unvisitedNum - 1;
            colony.ant(i).tour(j,vehTourLen(j)) = nextNode;

            [timeSlack, reservation] = resolveConflict(reservation, A, C, currentNode, nextNode, colony.ant(i).tick(j), implicitRoute, j);
            colony.ant(i).tick(j) = colony.ant(i).tick(j) + C(currentNode,nextNode) + timeSlack;

            colony.ant(i).tickHistory(j,vehTourLen(j)) = colony.ant(i).tick(j);

            % local pheromone update
            tau(currentNode,nextNode) = tau(currentNode,nextNode) * (1-psi) + tau0 * psi;
%         end
    end
    
    colony.ant(i).vehTourLen = vehTourLen;
    colony.ant(i).reservation = reservation;
    % complete the tour 
    for j = 1:vehNum
        colony.ant(i).tour(j,vehTourLen(j)+1) = homeIdx;
        colony.ant(i).tick(j,vehTourLen(j)+1) = ...
            colony.ant(i).tickHistory(j,vehTourLen(j)) + C(colony.ant(i).tour(j,vehTourLen(j)),homeIdx);
    end
end

end