function [ fitness, fitnessL, fitnessPerV, fitnessM, fitnessDev ] = fitnessFunctionVRPDM (tour, tourLen, graph, vehNum, implicitRoute)
global fitnessRes costWeight
N = graph.n;
nodePos = graph.node;
C = graph.edges;

fitnessPerVeh = zeros(vehNum,1);

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

fitnessM = max(fitnessPerVeh);
fitnessL = sum(fitnessPerVeh);
fitnessPerV = fitnessPerVeh;

posInfo = zeros(vehNum,fitnessRes,3);
devInfo = zeros(1,fitnessRes);
tStamp = linspace(1,fitnessM,fitnessRes);
for v = 1:vehNum
    detailTour = routeDecoder(nonzeros(tour(v,:)),implicitRoute, N);
    for t = 1:fitnessRes
        time = tStamp(t);
        posInfo(v,t,:) = getVehPos(detailTour, time, C, nodePos, implicitRoute, N);
    end
end

for i = 1:fitnessRes
    devInfo(i) = sum(std(reshape(posInfo(:,i,:),[vehNum,3,1])));
end

fitnessDev = sum(devInfo);

fitness = fitnessL - fitnessDev * costWeight;

end