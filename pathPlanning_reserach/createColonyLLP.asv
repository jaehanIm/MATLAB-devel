function [ colony ] = createColonyLLP( graph, colony, antNo, tau_tot, eta_tot, alpha,  beta)

global debugTemp

nodeNo = graph.n; %total Number of nodes
homeIdx = 1;
vehNum = 1;
cluNum = size(graph.totSubProbNodeIdx,2);
subNodeNo = [];
for i = 1:cluNum
    subNodeNo(i) = size(graph.totSubProbNodeIdx{i},1);
end
colonyHistory = [];

for i = 1 : antNo
    % prevEndNode = ??? for every ant.
    colony.ant(i).tour(:,1) = ones(vehNum,1) * 1; % starts from first node
    for c = 1:cluNum
        % CAUTION : Local idx used. 
        unvisitedNum = subNodeNo(c);
        unvisitedNum = unvisitedNum - 1;
        vehStepLen = ones(vehNum,1);
        
        if c == 1
            prevEndNode = []; % else, it will be designated in the previous loop
        end

        currSubNodes = graph.totSubProbNodeIdx{c};
        currSubNodes = vertcat(prevEndNode,currSubNodes);
        currEndNodes = graph.subProbEndNodeIdx{c};
        [~,endIdx] = ismember(currEndNodes,currSubNodes);

        tau = tau_tot(currSubNodes,currSubNodes);
        eta = eta_tot(currSubNodes,currSubNodes);

        while unvisitedNum ~= 0
            for j = 1:vehNum
                if unvisitedNum == 0
                    break;
                end
                % i = antNo, j = vehNum
                currentNode = colony.ant(i).tour(j,vehStepLen(j));
                P_allNodes = tau(currentNode,:).^alpha.*eta(currentNode,:).^beta;
                P_allNodes(nonzeros(colony.ant(i).tour(:))) = 0;
                P = P_allNodes./sum(P_allNodes);
                
                % for debugging
                debugTemp.P_allNodes = P_allNodes;
                debugTemp.currentNode = currentNode;
                debugTemp.i = i;
                debugTemp.j = j;
                debugTemp.colony = colony;
                debugTemp.vehTourLen = vehStepLen;

                % choose next node
                nextNode = rouletteWheel(P);
                vehStepLen(j) = vehStepLen(j) +1;
                unvisitedNum = unvisitedNum - 1;
                
                colony.ant(i).tour(j,vehStepLen(j)) = nextNode;
            end
        end
        
        colony.ant(i).vehTourLen = vehStepLen;

        % complete the tour by linking end node
        for j = 1:vehNum
            currentNode = colony.ant(i).tour(j,vehStepLen(j));
            if ~ismember(currentNode,endIdx)
                tau_end = tau([currentNode;endIdx],[currentNode;endIdx]);
                eta_end = eta([currentNode;endIdx],[currentNode;endIdx]);
                P_allNodes = tau_end(1,:).^alpha.*eta_end(1,:).^beta;
                P_allNodes(P_allNodes == inf) = 0;
                P = P_allNodes./sum(P_allNodes);
                finalNode = rouletteWheel(P);
                colony.ant(i).tour(j,vehStepLen(j)+1) = endIdx(finalNode-1);
            end % end of tour edition
        end % end of vehicle 
    end % end of cluster
    
    % decipher and add
    if c == 1
        colonyHistory = colony;
        for k = 1:antNo
            colonyHistory.ant(k).tour = graph.totSubProbNodeIdx{c}(colony.ant(k).tour);
        end
    else
        colonyHistory.ant(k).tour = vertcat(colonyHistory.ant(k).tour,map.totSubProbNodeIdx{c}(colony.ant(i).tour));
    end

end % end of iteration

end