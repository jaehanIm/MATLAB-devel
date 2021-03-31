function out = detectSubtours(graph)

n = size(graph,1);
is_visited = zeros(n,1);

startIdx = find(graph(:,1));
startNum = length(startIdx); % number of edges departing home node

subtourNum = 0;
subtourPath = [];

% From Home subtour detection
for i = 1:startNum
    curNode = 1;
    nextNode = 0;
    while nextNode ~= 1
        is_visited(curNode) = 1;
        if curNode == 1
            nextNode = startIdx(i);
        else
            nextNode = find(graph(:,curNode));
        end
        curNode = nextNode;
    end
end

% Undesired subtour detection
while ~isempty(find(~is_visited,1))
    subtourNum = subtourNum + 1;
    nextNode = 0;
    startNode = find(~is_visited,1,'first');
    curNode = startNode;
    count = 1;
    while nextNode ~= startNode
        subtourPath(subtourNum,count) = curNode;
        is_visited(curNode) = 1;
        nextNode = find(graph(:,curNode));
        curNode = nextNode;
        count = count + 1;
    end
end

out = subtourPath;

end