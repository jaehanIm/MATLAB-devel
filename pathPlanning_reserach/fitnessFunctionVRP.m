function [ fitness ] = fitnessFunctionVRP (tour, tourLen, graph, vehNum)

fitnessPerVeh = zeros(vehNum,1);

for i = 1:vehNum
    for j = 1:tourLen(i)
        currentNode = tour(i,j);
        nextNode = tour(i,j+1);
        fitnessPerVeh(i) = fitnessPerVeh(i) + graph.edges(currentNode, nextNode);
    end
end

fitness = max(fitnessPerVeh);

end