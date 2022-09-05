function [ colony ] = createColonyVRPTW( graph, colony, antNo, tau, eta, alpha,  beta, gamma, lambda, q, psi, vehNum, timeWindow, capacity, servTime)

global homeIdx debugTemp tau0

nodeNo = graph.n;
C = graph.edges;

for i = 1 : antNo
    
    unvisitedNum = nodeNo;
    
    initial_node = homeIdx;
    colony.ant(i).tour(:,1) = ones(vehNum,1) * initial_node;
    colony.ant(i).tick = zeros(vehNum,1);
    colony.ant(i).tickHistory = zeros(vehNum,1);
    unvisitedNum = unvisitedNum - 1;
    vehTourLen = ones(vehNum,1);
    colony.ant(i).tick = zeros(vehNum,1);

    while unvisitedNum ~= 0
%         for j = 1:vehNum
%             j = leastTourAgent(colony.ant(i).tour,C,vehNum,vehTourLen);
            j = getActiveVeh(colony.ant(i).tick, capacity);

            if unvisitedNum == 0
                break;
            end
            
            % i = antNo, j = vehNum
            currentNode = colony.ant(i).tour(j,vehTourLen(j));
            theta = updateTheta(nodeNo, C, currentNode, timeWindow, colony.ant(i).tick(j), lambda);

            P_allNodes = tau(currentNode,:).^alpha.*eta(currentNode,:).^beta.*theta.^gamma;
            P_allNodes(nonzeros(colony.ant(i).tour(:))) = 0;
            P = P_allNodes./sum(P_allNodes);
            
            % for debugging
            debugTemp.P_allNodes = P_allNodes;
            debugTemp.currentNode = currentNode;
            debugTemp.i = i;
            debugTemp.j = j;
            debugTemp.colony = colony;
            debugTemp.vehTourLen = vehTourLen;
            
            if rand(1) > q
                nextNode = rouletteWheel(P);
            else
                [~,nextNode] = max(P_allNodes);
                nextNode = nextNode(1); % just in case multiple nextNodes are generated
            end

            vehTourLen(j) = vehTourLen(j) +1;
            unvisitedNum = unvisitedNum - 1;
            colony.ant(i).tour(j,vehTourLen(j)) = nextNode;

            if timeWindow(nextNode,1) > colony.ant(i).tick(j) + C(currentNode,nextNode)
                colony.ant(i).tick(j) = timeWindow(nextNode,1) + servTime(nextNode);
            else
                colony.ant(i).tick(j) = colony.ant(i).tick(j) + C(currentNode,nextNode) + servTime(nextNode);
            end

            colony.ant(i).tickHistory(j,vehTourLen(j)) = colony.ant(i).tick(j);

            % local pheromone update
            tau(currentNode,nextNode) = tau(currentNode,nextNode) * (1-psi) + tau0 * psi;
%         end
    end
    
    colony.ant(i).vehTourLen = vehTourLen;
    
    % complete the tour 
    for j = 1:vehNum
        colony.ant(i).tour(j,vehTourLen(j)+1) = homeIdx;
    end
end

end