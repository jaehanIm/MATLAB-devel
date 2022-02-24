function [ fitness ] = fitnessFunctionLLP (tour, tourLen, graph)

fitness = 0;
vehNum = 1;

for i = 1:tourLen-1
    currentNode = tour(i);
    nextNode = tour(i+1);
    fitness = fitness + graph.edges(currentNode, nextNode);
end

end