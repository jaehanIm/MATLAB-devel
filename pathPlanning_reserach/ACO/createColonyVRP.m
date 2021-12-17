function [ colony ] = createColonyVRP( graph, colony, antNo, tau, eta, alpha,  beta)

global homeIdx vehNum temp mutationRate

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
            
            % for debugging
            temp.P_allNodes = P_allNodes;
            temp.currentNode = currentNode;
            temp.i = i;
            temp.j = j;
            temp.colony = colony;
            temp.vehTourLen = vehTourLen;
            
            if rand(1) < mutationRate
               P(P~=0) = 1;
            end
            
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