function [ colony ] = createColonyLLP( graph, colony, antNo, tau, eta, alpha,  beta)

global debugTemp

nodeNo = graph.n;
homeIdx = graph.homeIdx;
vehNum = graph.vehNum;

for i = 1 : antNo
    
    unvisitedNum = nodeNo;
    
    initial_node = homeIdx;
    colony.ant(i).tour(:,1) = ones(vehNum,1) * initial_node;
    unvisitedNum = unvisitedNum - 1;
    vehTourLen = ones(vehNum,1);

    while unvisitedNum ~= 0
        for j = 1:vehNum
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
            
            nextNode = rouletteWheel(P);
            vehTourLen(j) = vehTourLen(j) +1;
            unvisitedNum = unvisitedNum - 1;
            
            colony.ant(i).tour(j,vehTourLen(j)) = nextNode;
        end
    end
    
    colony.ant(i).vehTourLen = vehTourLen;
    
    % complete the tour by linking end node
    for j = 1:vehNum
        colony.ant(i).tour(j,vehTourLen(j)+1) = homeIdx;
    end
end

end