function colony = createColonyVRPMACC( graph, colony, antNo, tau, eta, param, vehNum)

global homeIdx debugTemp tau0

alpha = param.alpha;
beta = param.beta;
q = param.q;
psi = param.psi;

nodeNo = graph.n;
L = graph.edges;

for i = 1 : antNo
    
    unvisitedNum = nodeNo;
    
    initial_node = homeIdx;
    colony.ant(i).tour(:,1) = ones(vehNum,1) * initial_node;
    unvisitedNum = unvisitedNum - 1;
    vehTourLen = ones(vehNum,1);

    while unvisitedNum ~= 0
%         for j = 1:vehNum
            % tour construction / % i = antNo, j = vehNum
            j = leastTourAgent(colony.ant(i).tour,L,vehNum,vehTourLen);
            if unvisitedNum == 0
                break;
            end 
            
            % probability update
            currentNode = colony.ant(i).tour(j,vehTourLen(j));
            P_allNodes = tau(currentNode,:).^alpha.*eta(currentNode,:).^beta;
            P_allNodes(nonzeros(colony.ant(i).tour(:))) = 0;
            P = P_allNodes./sum(P_allNodes);
            
            % node addition / exploration vs exploitation
            if rand(1) > q
                nextNode = rouletteWheel(P); % exploration
            else
                [~,nextNode] = max(P_allNodes); % exploitation
                nextNode = nextNode(1); % just in case multiple nextNodes are generated
            end

            % information update
            vehTourLen(j) = vehTourLen(j) +1;
            unvisitedNum = unvisitedNum - 1;
            colony.ant(i).tour(j,vehTourLen(j)) = nextNode;

            % local pheromone update
            tau(currentNode,nextNode) = tau(currentNode,nextNode) * (1-psi) + tau0 * psi;

            % for debugging
            debugTemp.P_allNodes = P_allNodes;
            debugTemp.currentNode = currentNode;
            debugTemp.i = i;
            debugTemp.j = j;
            debugTemp.colony = colony;
            debugTemp.vehTourLen = vehTourLen;
%         end
    end
    
    colony.ant(i).vehTourLen = vehTourLen;
    
    % complete the tour 
    for j = 1:vehNum
        colony.ant(i).tour(j,vehTourLen(j)+1) = homeIdx;
    end
end

end