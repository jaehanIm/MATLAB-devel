function out = detectRoutes(graph)

n = size(graph,1);

startIdx = find(graph(:,1));
startNum = length(startIdx); % number of edges departing home node

subtourNum = 0;
subtourPath = [];

pathResult = [];
routeNum = 0;

% Route Detection
for i = 1:startNum
    curNode = 1;
    nextNode = 0;
    nodeNum = 1;
    while nextNode ~= 1
        if curNode == 1
            nextNode = startIdx(i);
        else
            nextNode = find(graph(:,curNode));
        end
        pathResult(i,nodeNum) = curNode;
        curNode = nextNode;
        nodeNum = nodeNum + 1;
    end
end

out = pathResult;

end