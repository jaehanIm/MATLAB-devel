function [ fitness, fitnessL, fitnessPerV, fitnessM, violation] = fitnessFunctionVRPTW (tour, tourLen, graph, vehNum, timeWindow, tickHistory, servTime, penalFactor)

fitnessPerVeh = zeros(vehNum,1);
violation = zeros(vehNum,1);

for i = 1:vehNum
    for j = 1:tourLen(i)
        currentNode = tour(i,j);
        nextNode = tour(i,j+1);
        fitnessPerVeh(i) = fitnessPerVeh(i) + graph.edges(currentNode, nextNode);
        if currentNode ~= 1
            initPenalty = max(timeWindow(currentNode,1) - (tickHistory(i,j) - servTime(currentNode)),0);
            endPenalty = max(tickHistory(i,j)-timeWindow(currentNode,2),0);
            violation(i) = violation(i) + initPenalty + endPenalty;
            fitnessPerVeh(i) = fitnessPerVeh(i) + (initPenalty + endPenalty)^2 * penalFactor;
        end
    end
end

fitnessM = max(fitnessPerVeh);
fitnessL = sum(fitnessPerVeh);
fitnessPerV = fitnessPerVeh;
fitness = fitnessL;

end