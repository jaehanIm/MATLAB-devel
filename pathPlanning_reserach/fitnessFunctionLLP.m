function [ fitness ] = fitnessFunctionLLP (tour, tourLen, graph)

fitness = 0;
vehNum = graph.vehNum;

for i = 1:vehNum
    for j = 1:tourLen(i)
        currentNode = tour(i,j);
        nextNode = tour(i,j+1);
        fitness = fitness + graph.edges(currentNode, nextNode);
    end
end

end