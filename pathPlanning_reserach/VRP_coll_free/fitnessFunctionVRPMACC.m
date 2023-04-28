function [ fitness, fitnessL, fitnessPerV, fitnessM ] = fitnessFunctionVRPMACC (tour, tourLen, graph, vehNum)

baseWeight = 1;
fitnessPerVeh = zeros(vehNum,1);

for i = 1:vehNum
    for j = 1:tourLen(i)
        if tourLen(i) ~= 1
            currentNode = tour(i,j);
            nextNode = tour(i,j+1);
            penaltyFactor = (baseWeight + j);
            penaltyFactor = 1;
            fitnessPerVeh(i) = fitnessPerVeh(i) + graph.edges(currentNode, nextNode) * penaltyFactor;
        else
            fitnessPerVeh(i) = 0;
        end
    end
end

fitnessM = max(fitnessPerVeh);
fitnessL = sum(fitnessPerVeh);
fitnessPerV = fitnessPerVeh;
fitness = fitnessL;

end