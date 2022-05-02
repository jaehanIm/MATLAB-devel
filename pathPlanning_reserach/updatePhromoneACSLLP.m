function [ tau ] = updatePhromoneACSLLP(tau , colony, rho, minIdx)

% Update the pheromone matrix
j = 1;
for k = 1:colony.ant(minIdx).vehStepLen - 1
    currentNode = colony.ant(minIdx).tour(k);
    nextNode = colony.ant(minIdx).tour(k+1);
    
    tau(currentNode , nextNode) = tau(currentNode , nextNode) * (1-rho)  + 1./ colony.ant(minIdx).fitness * rho;
    tau(nextNode , currentNode) = tau(nextNode , currentNode) * (1-rho)  + 1./ colony.ant(minIdx).fitness * rho;
end

% TODO : modify ACS udpate rule

end