function [newA,newC]=graphSparseConnection(node,A,C,L)

while true
G = graph(C);
discon = conncomp(G);
disconNum = max(conncomp(G));
centerPos = [];
if disconNum > 1
%     disp("Disconnection Detected!")
    % CG
    for i = 1:disconNum
        centerPos(i,:) = mean(node(find(discon==i),:),1);
    end

    disconBridgeLen = sparse(disconNum,disconNum);
    for i = 1:disconNum
        for j = 1:disconNum
            disconBridgeLen(i,j) = norm(centerPos(i,:) - centerPos(j,:));
        end
        disconBridgeLen(i,i) = 1e10;
    end

    [~,I]=min(disconBridgeLen,[],'all');
    [group1,group2]=ind2sub(size(disconBridgeLen),I);

    group1Nodes = find(discon==group1);
    group2Nodes = find(discon==group2);

    L_interest = L(group1Nodes,group2Nodes);
    [~,I]=min(L_interest,[],'all');
    [row,col] = ind2sub(size(L_interest),I);
    
    A(group1Nodes(row),group2Nodes(col)) = 1;
    A(group2Nodes(col),group1Nodes(row)) = 1;

    C(group1Nodes(row),group2Nodes(col)) = L(group1Nodes(row),group2Nodes(col));
    C(group2Nodes(col),group1Nodes(row)) = L(group1Nodes(row),group2Nodes(col));

else
    break;
end
end

newA = A;
newC = C;

end