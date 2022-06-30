function [route,routeL]=NNHeuristic(map)

% require N(number of nodes), node, full length matrix L
N = map.N;
net = map.node;
costMatrix = map.L;
costMatrix(costMatrix == 0) = nan;

% home, end point initialization
record = zeros(N,size(net,2));
record(1,:) = net(1,:);
record(end,:) = net(end,:);

% run NN solver
cost_NNM = 0;
curridx = 1;
costMatrix_temp = costMatrix(:,1:end-1);
costMatrix_temp(:,1) = nan;
for i = 1:length(net)-2
[~,j]=min(costMatrix_temp(curridx,:));
record(i+1,:) = net(j,:);
cost_NNM = cost_NNM + costMatrix_temp(curridx,j);
curridx = j;
costMatrix_temp(:,j) = nan;
end
% cost_NNM = cost_NNM + costMatrix(curridx,end); %restricted
% record = vertcat(record,net(end,:)); %restricted

route = record;
routeL = cost_NNM;

if N == 1
    routeL = 0;
elseif N <= 2
    routeL = costMatrix(1,2);
end

end