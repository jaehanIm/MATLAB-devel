function [ fitness, fitnessL, fitnessPerV, fitnessM, violation] = fitnessFunctionVRPCF (tour, tourLen, graph, vehNum, tickHistory, servTime, param)

fitnessPerVeh = zeros(vehNum,1);
violation = zeros(vehNum,1);

for i = 1:vehNum
    for j = 1:tourLen(i)
        currentNode = tour(i,j);
        nextNode = tour(i,j+1);
        fitnessPerVeh(i) = fitnessPerVeh(i) + graph.edges(currentNode, nextNode);
    end
end

fitnessM = max(fitnessPerVeh);
fitnessL = sum(fitnessPerVeh);
fitnessPerV = fitnessPerVeh;
fitness = fitnessL;

end