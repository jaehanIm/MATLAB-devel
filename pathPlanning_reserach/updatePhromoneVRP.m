function [ tau ] = updatePhromoneVRP(tau , colony)
global vehNum

% Update the pheromone matrix
for i = 1:length(colony.ant)
    for j = 1:vehNum
        for k = 1:colony.ant(i).vehTourLen(j) - 1
            currentNode = colony.ant(i).tour(j,k);
            nextNode = colony.ant(i).tour(j,k+1);
            
            tau(currentNode , nextNode) = tau(currentNode , nextNode)  + 1./ colony.ant(i).fitness;
            tau(nextNode , currentNode) = tau(nextNode , currentNode)  + 1./ colony.ant(i).fitness;

        end
    end
end

end