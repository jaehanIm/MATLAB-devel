function [ colony ] = createColonyVRP( graph, colony, antNo, tau, eta, alpha,  beta)

global homeIdx vehNum temp temp2

nodeNo = graph.n;

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
            currentNode = colony.ant(i).tour(j,vehTourLen(j));
            P_allNodes = tau(currentNode,:).^alpha.*eta(currentNode,:).^beta;
            P_allNodes(nonzeros(colony.ant(i).tour(:))) = 0;
            P = P_allNodes./sum(P_allNodes);
            
            temp{1} = P_allNodes;
            temp{2} = currentNode;
            temp{3} = i;
            temp{4} = j;
            temp{5} = colony;
            temp{6} = vehTourLen;
            
            nextNode = rouletteWheel(P);
            vehTourLen(j) = vehTourLen(j) +1;
            unvisitedNum = unvisitedNum - 1;
            
            colony.ant(i).tour(j,vehTourLen(j)) = nextNode;
        end
    end
    
    colony.ant(i).vehTourLen = vehTourLen;
    
    % complete the tour 
    for j = 1:vehNum
        colony.ant(i).tour(j,vehTourLen(j)+1) = homeIdx;
    end
end

end