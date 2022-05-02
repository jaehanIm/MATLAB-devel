function [ tau ] = updatePhromoneACSVRP(tau , colony, rho, minIdx, vehNum)

% Update best ant pheromone trail
for j = 1:vehNum
    for k = 1:colony.ant(minIdx).vehTourLen(j) - 1
        currentNode = colony.ant(minIdx).tour(j,k);
        nextNode = colony.ant(minIdx).tour(j,k+1);
        
        tau(currentNode , nextNode) = tau(currentNode , nextNode) * (1-rho)  + 1./ colony.ant(minIdx).fitness * rho;
        tau(nextNode , currentNode) = tau(nextNode , currentNode) * (1-rho)  + 1./ colony.ant(minIdx).fitness * rho;
    end
end

end