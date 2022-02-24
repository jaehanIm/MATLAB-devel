function [ tau ] = updatePhromoneLLP(tau , colony)

% Update the pheromone matrix
for i = 1:length(colony.ant)
    for k = 1:colony.ant(i).vehStepLen - 1
        currentNode = colony.ant(i).tour(k);
        nextNode = colony.ant(i).tour(k+1);
        
        tau(currentNode , nextNode) = tau(currentNode , nextNode)  + 1./ colony.ant(i).fitness;
        tau(nextNode , currentNode) = tau(nextNode , currentNode)  + 1./ colony.ant(i).fitness;
    end
end

% TODO : modify ACS udpate rule

end