function [ colonyHistory ] = createColonyACSLLP( graph, colony, antNo, tau_tot, eta_tot, alpha,  beta, q, psi)

global debugTemp tau0

nodeNo = graph.n; %total Number of nodes
homeIdx = 1;
vehNum = 1;
cluNum = size(graph.totSubProbNodeIdx,2);
subNodeNo = [];
for i = 1:cluNum
    subNodeNo(i) = size(graph.totSubProbNodeIdx{i},1);
end
colonyHistory = [];
colonyCluster = [];

bridgeCost = zeros(cluNum,1);

for i = 1 : antNo
    for c = 1:cluNum
        % CAUTION : Local idx used. 
        colony = [];
        unvisitedNum = subNodeNo(c);
        
        if c == 1
            prevEndNode = []; % else, it will be designated in the previous loop
            unvisitedNum = unvisitedNum - 1;
        end

        colony.tour(:,1) = ones(vehNum,1) * 1; % starts from first node
        vehStepLen = ones(vehNum,1);

        currSubNodes = graph.totSubProbNodeIdx{c};
        currSubNodes = vertcat(prevEndNode,currSubNodes);
        currEndNodes = graph.subProbEndNodeIdx{c};
        [~,endIdx] = ismember(currEndNodes,currSubNodes);

        tau = tau_tot(currSubNodes,currSubNodes);
        eta = eta_tot(currSubNodes,currSubNodes);
        C = graph.edges(currSubNodes,currSubNodes);

        while unvisitedNum ~= 0
            for j = 1:vehNum
                if unvisitedNum == 0
                    break;
                end
                % i = antNo, j = vehNum
                currentNode = colony.tour(j,vehStepLen(j));
                P_allNodes = tau(currentNode,:).^alpha.*eta(currentNode,:).^beta;
                P_allNodes(nonzeros(colony.tour(:))) = 0;
                P = P_allNodes./sum(P_allNodes);
                
                % for debugging
                debugTemp.P_allNodes = P_allNodes;
                debugTemp.currentNode = currentNode;
                debugTemp.i = i;
                debugTemp.j = j;
                debugTemp.colony = colony;
                debugTemp.vehTourLen = vehStepLen;

                % choose next node
                if rand(1) > q
                    nextNode = rouletteWheel(P);
                else
                    [~,nextNode] = max(P_allNodes);
                    nextNode = nextNode(1); % just in case multiple nextNodes are generated
                end

                vehStepLen(j) = vehStepLen(j) +1;
                unvisitedNum = unvisitedNum - 1;
                
                colony.tour(j,vehStepLen(j)) = nextNode;

                % local pheromone update
                tau(currentNode,nextNode) = tau(currentNode,nextNode) * (1-psi) + tau0 * psi;
            end
        end

        % complete the tour by linking end node
        for j = 1:vehNum
            currentNode = colony.tour(j,vehStepLen(j));
            if ~ismember(currentNode,endIdx)
                tau_end = tau([currentNode;endIdx],[currentNode;endIdx]);
                eta_end = eta([currentNode;endIdx],[currentNode;endIdx]);
                P_allNodes = tau_end(1,:).^alpha.*eta_end(1,:).^beta;
                P_allNodes(P_allNodes == inf) = 0;
                P = P_allNodes./sum(P_allNodes);
                finalNode = rouletteWheel(P);
                colony.tour(j,vehStepLen(j)+1) = endIdx(finalNode-1);
                vehStepLen(j) = vehStepLen(j) + 1;
            end % end of tour edition
        end % end of vehicle 
        colony.vehStepLen = vehStepLen;

        % decipher and add to history
        if c == 1
            colonyCluster.tour = currSubNodes(colony.tour);
            colonyCluster.vehStepLen = colony.vehStepLen;
        else
            colonyCluster.tour = vertcat(colonyCluster.tour, currSubNodes(colony.tour));
            colonyCluster.vehStepLen = colonyCluster.vehStepLen + colony.vehStepLen;
        end

        prevEndNode = currSubNodes(colony.tour(end));
        bridgeCost(c) = graph.edges(currSubNodes(colony.tour(1)),currSubNodes(colony.tour(2)));
        clusterCost(c) = fitnessFunctionLLP(currSubNodes(colony.tour(2:end)),colony.vehStepLen-1,graph);

    end % end of cluster
    
    % add depot node
    colonyCluster.tour(end+1) = 1;
    colonyCluster.vehStepLen = colonyCluster.vehStepLen + 1;
    colonyCluster.bridgeCost = bridgeCost;
    colonyCluster.clusterCost = clusterCost;
    colonyHistory.ant(i) = colonyCluster;
end % end of all ants

end