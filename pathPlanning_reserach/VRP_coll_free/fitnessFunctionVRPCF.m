function [ fitness, fitnessL, fitnessPerV, fitnessM, violation] = fitnessFunctionVRPCF (tour, tourLen, graph, vehNum, tickHistory, servTime, param, incompleteFlag, N)


fitnessPerVeh = zeros(vehNum,1);
violation = zeros(vehNum,1);

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

fitnessE = zeros(vehNum,1);
for i = 1:vehNum
    fitnessE(i) = tickHistory(i,tourLen(i));
end
fitnessE = max(fitnessE);
fitnessM = max(fitnessPerVeh);
fitnessL = sum(fitnessPerVeh);
fitnessPerV = fitnessPerVeh;

fitness = fitnessM;

coveredNum = length(unique(tour(:)));
if coveredNum ~= N+1
    disp("Infeasible! cursed!")
    fitness = fitness * exp(N+1-coveredNum);
end

end