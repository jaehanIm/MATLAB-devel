function colony = createColonyVRPDM( graph, colony, antNo, tau, eta, param, vehNum, implicitRoute, nodePos)

global homeIdx debugTemp tau0

alpha = param.alpha;
beta = param.beta;
gamma = param.gamma;
q = param.q;
psi = param.psi;

N = graph.n;
L = graph.edges; %NxN cost matrix

for i = 1 : antNo
    
    unvisitedNum = N;
    
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
            P_allNodes = P_allNodes./sum(P_allNodes);

%             P_dev = devMaxPUpdate(colony.ant(i).tour, implicitRoute, currentNode, P_allNodes, nodePos, L, j, N);
%             P_dev = 1./P_dev;
%             P_dev(find(isinf(P_dev))) = 0;
%             P_dev = P_dev./sum(P_dev);

            P_dev = 0;

            P = P_allNodes + P_dev *gamma ;
            P = P./sum(P);
            
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
            debugTemp.P_dev = P_dev;
            debugTemp.P = P;
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