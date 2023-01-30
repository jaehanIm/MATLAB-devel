function [ fitness, fitnessL, fitnessPerV, fitnessM ] = fitnessFunctionVRP (tour, tourLen, graph, vehNum)


fitnessPerVeh = zeros(vehNum,1);
fitnessPerVeh_e = zeros(vehNum,1);

for i = 1:vehNum
    for j = 1:tourLen(i)
        if tourLen(i) ~= 1
            currentNode = tour(i,j);
            nextNode = tour(i,j+1);
            fitnessPerVeh(i) = fitnessPerVeh(i) + graph.edges(currentNode, nextNode);
        else
            fitnessPerVeh(i) = 0;
        end
    end
end

for i = 1:vehNum
    for j = 1:tourLen(i)-1
        currentNode = tour(i,j);
        nextNode = tour(i,j+1);
        fitnessPerVeh_e(i) = fitnessPerVeh_e(i) + graph.edges(currentNode, nextNode);
    end
end

fitnessM = max(fitnessPerVeh);
fitnessL = sum(fitnessPerVeh);
fitnessPerV = fitnessPerVeh;
fitnessME = max(fitnessPerVeh_e);
fitness = fitnessM;

end