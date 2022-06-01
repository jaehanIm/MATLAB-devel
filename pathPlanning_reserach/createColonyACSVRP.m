function [ colony ] = createColonyACSVRP( graph, colony, antNo, tau, eta, alpha,  beta, q, psi, vehNum)

global homeIdx debugTemp tau0

nodeNo = graph.n;
L = graph.edges;

for i = 1 : antNo
    
    unvisitedNum = nodeNo;
    
    initial_node = homeIdx;
    colony.ant(i).tour(:,1) = ones(vehNum,1) * initial_node;
    unvisitedNum = unvisitedNum - 1;
    vehTourLen = ones(vehNum,1);

    while unvisitedNum ~= 0
        for j = 1:vehNum
%             j = leastTourAgent(colony.ant(i).tour,L,vehNum,vehTourLen);
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

            % local pheromone update
            tau(currentNode,nextNode) = tau(currentNode,nextNode) * (1-psi) + tau0 * psi;
        end
    end
    
    colony.ant(i).vehTourLen = vehTourLen;
    
    % complete the tour 
    for j = 1:vehNum
        colony.ant(i).tour(j,vehTourLen(j)+1) = homeIdx;
    end
end

end