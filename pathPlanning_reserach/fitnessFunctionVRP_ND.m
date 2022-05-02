function [ fitness ] = fitnessFunctionVRP_ND (tour, tourLen, graph, vehNum, ND)

fitnessPerVeh = zeros(vehNum,1);

for i = 1:vehNum
    for j = 1:tourLen(i)
        currentNode = tour(i,j);
        nextNode = tour(i,j+1);
        fitnessPerVeh(i) = fitnessPerVeh(i) + graph.edges(currentNode, nextNode) + ND(nextNode);
    end
end

fitness = max(fitnessPerVeh);

end