function P_dev = devMaxPUpdate(tour, implicitRoute, curNode, P_allNodes, nodePos, C, egoVIdx, N)
validNodes = find(P_allNodes~=0);
P_dev = zeros(1,N);

for i = 1:length(validNodes)
    nextNode = validNodes(i);
    P_dev(nextNode) = calcDev(tour, implicitRoute, curNode, nextNode, nodePos, C, egoVIdx, N);
end

end