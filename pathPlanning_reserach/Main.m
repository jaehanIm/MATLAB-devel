clear all

addpath('./ACO')
addpath('./MILP')
addpath('./ComDetTBv090/')
addpath('./ComDetTBv090/Algorithms/')
addpath('./ComDetTBv090/Auxiliary/')
addpath('./ComDetTBv090/Graphs/')

%% Param setting
vnum = 3;
depotPos = [60 60 10];
fovFactor = 3;
inpection_dist = 7; % Inspection distance
mapheight = 4.0;
conThres = 20;
antNo = 20;
stopThres = 200;

%% wp generator
poc_path_planner

node = [airPosX(~isnan(airPosZ(:))),airPosY(~isnan(airPosZ(:))),airPosZ(~isnan(airPosZ(:)))];
node = vertcat(depotPos,node); % add depot

N = size(node,1);
A = zeros(N,N); % connectivity matrix
C = sparse(N,N); % cost matrix
L = sparse(N,N); % linear distance matrix

for i = 1:N
    for j = 1:N
        L(i,j) = norm(node(i,:)-node(j,:));
    end
end


%% Build Network Graph
for i = 2:N % do not connect depot!
    for j = 2:N
        if i~=j
            if L(i,j) < conThres
                A(i,j) = 1;
                C(i,j) = L(i,j);
            end
        else
            A(i,j) = 0;
        end
    end
end
[A,C]=graphSparseConnection(node,A,C,L);
A(1,:) = 0; A(:,1) = 0;

G = graph(C);
A_orig = A;
C_orig = C;

totalDegree = sum(A(2:end,2:end),'all')/2;
completeDegree = nchoosek(N-1,2);
degreeConnectivity = totalDegree / completeDegree


%% Save Network
graph1.n = N;
graph1.node.x = node(:,1);
graph1.node.y = node(:,2);
graph1.node.z = node(:,3);
graph1.edges = L;

save('graph_complete.mat','graph1');

%% Graph Clustering
tic
A = A_orig;
A_temp = A(2:end,2:end); %except home node
C_temp = C(2:end,2:end);


cluIdx = GCModulMax2(A_temp); % DO NOT USE GCModulMax1 -> erroneous algorithm
% cluIdx = GCSpectralClust1(C_temp,vnum*5); cluIdx = cluIdx(:,end);
% cluIdx = GCModulMax3(C_temp);
% cluIdx = spectralcluster(A_temp,vnum*2,'LaplacianNormalization','symmetric');

cluIdx = vertcat(1,cluIdx+1);
cluNum = max(cluIdx);

nodeInCluIdx = [];
for i = 1:cluNum
    nodeInCluIdx{i} = find(cluIdx == i);
end
intraCluIdxSet = [];
intraCluNodeSet = [];
for i = 1:cluNum-1
    for j = i+1:cluNum
        [I,J]=find(A(nodeInCluIdx{i},nodeInCluIdx{j}));
        intraCluIdxSet{i,j} = [I,J];
        intraCluIdxSet{j,i} = [J,I];
        intraCluNodeSet{i,j} = [nodeInCluIdx{i}(I),nodeInCluIdx{j}(J)];
        intraCluNodeSet{j,i} = [nodeInCluIdx{j}(J),nodeInCluIdx{i}(I)];
    end
end

% connect depot to cluster
for j = 2:cluNum
    [~,I]=min(L(1,nodeInCluIdx{j}));
    intraCluIdxSet{1,j} = [1,I];
    intraCluIdxSet{j,1} = [I,1];
    intraCluNodeSet{1,j} = [1,nodeInCluIdx{j}(I)];
    intraCluNodeSet{j,1} = [nodeInCluIdx{j}(I),1];
    C(1,nodeInCluIdx{j}(I)) = L(1,nodeInCluIdx{j}(I));
    C(nodeInCluIdx{j}(I),1) = L(nodeInCluIdx{j}(I),1);
    A(1,nodeInCluIdx{j}(I)) = 1;
    A(nodeInCluIdx{j}(I),1) = 1;
end

clusteringTime = toc;

figure(5)
clf
p = plot(graph(A),'Layout','force','EdgeAlpha',0.3,'MarkerSize',3);
p.NodeCData = cluIdx;
colormap jet
colorbar
drawnow

figure(52)
clf
p = plot(graph(A),'Layout','force','EdgeAlpha',0.3,'MarkerSize',3);


%% SuperNet Construction
superNet = [];
superNet.A = sparse(cluNum,cluNum);
superNet.C = sparse(cluNum,cluNum);
superNet.ND = sparse(cluNum,1); % super cluster demand (node demand)
superNet.pos = [];

for i = 1:cluNum-1
    for j = i+1:cluNum
        if ~isempty(intraCluIdxSet{i,j})
            superNet.A(i,j) = 1;
            superNet.A(j,i) = 1;
            tempL = sparse(1,size(intraCluNodeSet{i,j},1));
            for k = 1:size(intraCluNodeSet{i,j},1)
                tempL(k) = C(intraCluNodeSet{i,j}(k,1),intraCluNodeSet{i,j}(k,2));
            end
            superNet.C(i,j) = min(tempL);
            superNet.C(j,i) = min(tempL);
        end
    end
end

for i = 2:cluNum
    superNet.C(1,i) = norm(depotPos - mean(node(nodeInCluIdx{i},:)));
    superNet.C(i,1) = superNet.C(1,i);
end

% for i = 2:cluNum
%     tempL = [];
%     for j = 1:size(nodeInCluIdx{i},1)
%         tempL(j) = norm(depotPos - node(nodeInCluIdx{i},:));
%     end
%     superNet.C(1,i) = min(tempL);
%     superNet.C(i,1) = superNet.C(1,i);
% end

superNet.ND(1) = 1;
for i = 2:cluNum
    superNet.ND(i) = size(nodeInCluIdx{i},1);
end

for i = 1:cluNum
    if size(nodeInCluIdx{i},1) == 1
        superNet.pos(i,1:3) = node(nodeInCluIdx{i},:);
    else
        superNet.pos(i,1:3) = mean(node(nodeInCluIdx{i},:));
    end
end

superNet.C_temp = superNet.C;


figure(4)
clf
grid on
hold on
for i = 1:size(G.Edges,1)
    startIdx = G.Edges.EndNodes(i,1);
    EndIdx = G.Edges.EndNodes(i,2);
    startPos = node(startIdx,:);
    EndPos = node(EndIdx,:);
    line([startPos(1) EndPos(1)],[startPos(2) EndPos(2)],[startPos(3) EndPos(3)]);
end
% for i = 1:cluNum
%     temp = node(nodeInCluIdx{i},:);
%     plot3(temp(:,1),temp(:,2),temp(:,3),'x','LineWidth',5,'MarkerSize',5);
% end
plot3(node(:,1),node(:,2),node(:,3),'kx','LineWidth',5,'MarkerSize',5)
mesh(voxelPosX,voxelPosY,voxelFilterData);
% plot3(superNet.pos(:,1),superNet.pos(:,2),superNet.pos(:,3),'yx','MarkerSize',10,'LineWidth',5)
% for i = 1:cluNum
%     text(superNet.pos(i,1),superNet.pos(i,2),superNet.pos(i,3),num2str(i));
% end
% for i = 1:N
%     text(node(i,1),node(i,2),node(i,3),num2str(i));
% end

%% SuperNet (HLP) initialization
for i = 1:cluNum
    localNet.N = size(nodeInCluIdx{i},1);
    localNet.node = node(nodeInCluIdx{i},:);
    localNet.L = L(nodeInCluIdx{i},nodeInCluIdx{i});
    [~,localCost] = NNHeuristic(localNet);
    superNet.ND(i) = localCost;
end
superNet.ND(1) = 1; % depot is nan

%% Local Network (LLP) Completefication

% superNet complefication (initialization)
superPosC = sparse(cluNum,cluNum);
for i = 1:cluNum-1
    for j = i+1:cluNum
        superPosC(i,j) = norm(superNet.pos(i,:)-superNet.pos(j,:));
        superPosC(j,i) = superPosC(i,j);
    end
end

tic
disconCluList = []; disconCount = 1;
for i = 1:cluNum-1
    for j = i+1:cluNum
        if superNet.A(i,j) == 0
            disconCluList(disconCount,:) = [i,j];
            disconCount = disconCount + 1;
            disp("[Warning] disconnected cluster detected");  
            distInterest = L(nodeInCluIdx{i},nodeInCluIdx{j});
            [val,I] = min(distInterest,[],'all');
            [idxx,idxy] = ind2sub(size(distInterest),I);
            [implRoute_temp,d] = shortestpath(G,nodeInCluIdx{i}(idxx),nodeInCluIdx{j}(idxy),'Method','positive');
            superNet.C(i,j) = d;
            superNet.C(j,i) = d;
            implicitRoute{nodeInCluIdx{i}(idxx),nodeInCluIdx{j}(idxy)} = implRoute_temp;
            implicitRoute{nodeInCluIdx{j}(idxy),nodeInCluIdx{i}(idxx)} = fliplr(implRoute_temp);
            intraCluNodeSet{i,j} = [nodeInCluIdx{i}(idxx),nodeInCluIdx{j}(idxy)];
            intraCluNodeSet{j,i} = [nodeInCluIdx{j}(idxy),nodeInCluIdx{i}(idxx)];
            C(nodeInCluIdx{i}(idxx),nodeInCluIdx{j}(idxy)) = d;
            C(nodeInCluIdx{j}(idxy),nodeInCluIdx{i}(idxx)) = d;
        end
    end
end

interCompleteTime = toc;

tic
% inter-intra cluster complefication
implicitRoute = [];
for i = 1:cluNum-1
    localNodeIdx = nodeInCluIdx{i};
    localIntraNodeIdx = [];
    for j = i+1:cluNum
        if ~isempty(intraCluNodeSet{i,j})
            localIntraNodeIdx = intraCluNodeSet{i,j}(:,2);
            localIntraNodeIdx = unique(localIntraNodeIdx);
            locTotNodIdx = vertcat(localNodeIdx,localIntraNodeIdx);
            % start connection
            localC = C(locTotNodIdx,locTotNodIdx);
            localG = graph(localC);
            for ki = 1:size(locTotNodIdx,1)-1
                for kj = ki+1:size(locTotNodIdx,1)
                    origIdxi = locTotNodIdx(ki);
                    origIdxj = locTotNodIdx(kj);
                    if A(origIdxi,origIdxj) == 0
                        [implRoute,d]=shortestpath(localG,ki,kj,'Method','positive');
                        implicitRoute{origIdxi,origIdxj} = locTotNodIdx(implRoute);
                        implicitRoute{origIdxj,origIdxi} = fliplr(locTotNodIdx(implRoute));
                        if d == inf
                            disp("[CRITICAL] FAULTY CLUSTERING ALGORITHM - STOP IMMEDIATELY")
                        end
                        C(origIdxi,origIdxj) = d;
                        C(origIdxj,origIdxi) = d;
                    end
                end
            end
            % connection complete
        end
    end
end

implicitRoute = [];
for i = fliplr(2:cluNum)
    localNodeIdx = nodeInCluIdx{i};
    localIntraNodeIdx = [];
    for j = fliplr(1:i-1)
        if ~isempty(intraCluNodeSet{i,j})
            localIntraNodeIdx = intraCluNodeSet{i,j}(:,2);
            localIntraNodeIdx = unique(localIntraNodeIdx);
            locTotNodIdx = vertcat(localNodeIdx,localIntraNodeIdx);
            % start connection
            localC = C(locTotNodIdx,locTotNodIdx);
            localG = graph(localC);
            for ki = 1:size(locTotNodIdx,1)-1
                for kj = ki+1:size(locTotNodIdx,1)
                    origIdxi = locTotNodIdx(ki);
                    origIdxj = locTotNodIdx(kj);

                    if A(origIdxi,origIdxj) == 0
                        [implRoute,d]=shortestpath(localG,ki,kj,'Method','positive');
                        implicitRoute{origIdxi,origIdxj} = locTotNodIdx(implRoute);
                        implicitRoute{origIdxj,origIdxi} = fliplr(locTotNodIdx(implRoute));
                        C(origIdxi,origIdxj) = d;
                        C(origIdxj,origIdxi) = d;
                    end
                end
            end
            % connection complete
        end
    end
end

intraCompleteTime = toc;

figure(44)
clf
grid on
hold on
% for i = 1:size(G.Edges,1)
%     startIdx = G.Edges.EndNodes(i,1);
%     EndIdx = G.Edges.EndNodes(i,2);
%     startPos = node(startIdx,:);
%     EndPos = node(EndIdx,:);
%     line([startPos(1) EndPos(1)],[startPos(2) EndPos(2)],[startPos(3) EndPos(3)]);
% end
for i = 1:cluNum
    temp = node(nodeInCluIdx{i},:);
    plot3(temp(:,1),temp(:,2),temp(:,3),'x','LineWidth',5,'MarkerSize',5);
end
% mesh(voxelPosX,voxelPosY,voxelFilterData);
plot3(superNet.pos(:,1),superNet.pos(:,2),superNet.pos(:,3),'yx','MarkerSize',10,'LineWidth',5)
for i = 1:cluNum
    text(superNet.pos(i,1),superNet.pos(i,2),superNet.pos(i,3),num2str(i));
end
for i = 1:N
    text(node(i,1),node(i,2),node(i,3),num2str(i));
end
axis equal
drawnow

%% Solver
tic
firstTry = true;
prevSuperRoute = [];
superRoute = sparse(vnum,cluNum-vnum+1);
trialNum = 1;
totalScoreHistory = [];
totalScoreHistoryL = [];
totalTourHistory = [];
terminationType = [];
superRouteHistory = [];
vehScoreHistory = [];
while true
    disp(["==== Trial Num : "+num2str(trialNum)+" ===="]);

    if vnum >= cluNum
        disp("[error] cannot solve VRP : too many vehicles!")
        break;
    end
    
    %% HLP
%     HLP_complete = false;
%     map_H.A = superNet.A;
%     map_H.C = superNet.C;
%     map_H.ND = superNet.ND;
%     map_H.vnum = vnum;
%     map_H.N = cluNum; 
%     if ~firstTry
%         map_H.totLoad = totalScore/vnum*0.2;
%     else
%         map_H.totLoad = mean(superNet.C(2:end,2:end),'all')*(cluNum/vnum-1) + mean(superNet.ND(2:end))*cluNum/vnum;
%         map_H.totLoad = map_H.totLoad * 0.5;
%     end
% 
%     while ~HLP_complete
%         [superRoute,superScore]=HLP_solver(map_H);
%         if superScore == -1
%             map_H.totLoad = map_H.totLoad * 1.1;
%         else
%             HLP_complete = true;
%             break;
%         end
%     end

%     [superRoute,superScore]=HLP_solver(map_H);
%     map_H = [];
%     map_H.N = N;
%     map_H.C = C;
%     map_H.totSubProbNodeIdx = totSubProbNodeIdx;
%     map_H.subProbEndNodeIdx = subProbEndNodeIdx;
%     [tour] = ACOVRP_forSoleACO(map,antNo,termCond)
%     [tour,score,clusterCost,bridgeCost,residueCost]=LLP_solver_ACS(map,20,400);
    map_H.n = cluNum;
    map_H.edges = superNet.C;
    map_H.ND = superNet.ND;
    map_H.vnum = vnum;
    superRoute = HLP_solver_ACS(map_H, 20, 400);

    superRoute

    %% Check Termination Condition
    if ~firstTry
        %Check termination condition
        HLPSolIdentCheck = 0;
        if size(superRoute) == size(prevSuperRoute)
            curRouteGroup = sparse(vnum,size(superRoute,2));
%             curRouteGroup = [];
            for v = 1:vnum
                uniqueSuperRouteLen = length(unique(superRoute(v,:)));
                curRouteGroup(v,1:uniqueSuperRouteLen) = unique(superRoute(v,:));
            end
            for v = 1:vnum
                for vj = 1:vnum
                    prevRouteGroup = unique(prevSuperRoute(vj,:));
                    if ismember(curRouteGroup(v,:),prevRouteGroup)
                        HLPSolIdentCheck = HLPSolIdentCheck+1;
                        break;
                    end
                end
            end
            if HLPSolIdentCheck == vnum
                disp("Identical HLP solution derived.");
                disp("Termination condition met. Finishing the solver!");
                terminationType = "HLPCONV";
                break;
            end
        end
    end

    %% LLP
    scoreRecord = [];
    tourRecord = [];
    for v = 1:vnum % for each vehicle (super route set)
        totSubProbNodeIdx = [];
        subProbEndNodeIdx = [];
        superRouteL = sum((superRoute(v,:)~=0));
        for c = 2:superRouteL %for number of involved clusters            
            currClus = superRoute(v,c);
%             prevClus = superRoute(v,c-1);
            if c ~= superRouteL
                nextClus = superRoute(v,c+1);
            else
                nextClus = superRoute(v,1);
            end
            totSubProbNodeIdx{c-1} = nodeInCluIdx{currClus};
            if c == 2
                prevClusEndNodeIdx = 1;
                totSubProbNodeIdx{c-1} = vertcat(prevClusEndNodeIdx,totSubProbNodeIdx{c-1});
            end
            subProbEndNodeIdx{c-1} = intraCluNodeSet{currClus,nextClus}(:,1);
            subProbEndNodeIdx{c-1} = unique(subProbEndNodeIdx{c-1});
        end
%         totSubProbNodeIdx
%         subProbEndNodeIdx

        % fomulate data structure for ACO
        map = [];
        map.N = N;
        map.C = C;
        map.totSubProbNodeIdx = totSubProbNodeIdx;
        map.subProbEndNodeIdx = subProbEndNodeIdx;
            
        % run ACO
%         [tour,score,clusterCost,bridgeCost,residueCost]=LLP_solver(map,antNo,stopThres);
        [tour,score,clusterCost,bridgeCost,residueCost]=LLP_solver_ACS(map,antNo,stopThres);
        scoreRecord(v) = score;       
        tourRecord{v} = tour;
        costRecord(v).clusterCost = clusterCost;
        costRecord(v).bridgeCost = bridgeCost;
        costRecord(v).residueCost = residueCost;

        % update super network
    end
    totalScoreL = sum(scoreRecord);
    totalScore = max(scoreRecord);
    vehScore = scoreRecord;
    solveTime(trialNum) = toc;

    %% Update Super Network
    for v = 1:vnum
        superRouteL = sum((superRoute(v,:)~=0));
        
        % initRoute update
        initClu = superRoute(v,2);
        superNet.C(1,initClu) = costRecord(v).residueCost(1); % initCost
        superNet.C(initClu,1) = superNet.C(1,initClu);

        % cluster ND update
        for i = 1:superRouteL-1
            superNet.ND(superRoute(v,i+1)) = costRecord(v).clusterCost(i);
        end

        % bridge update
        if ~isempty(costRecord(v).bridgeCost)
            for i = 1:superRouteL-2
                firstClu = superRoute(v,i+1);
                secondClu = superRoute(v,i+2);
                superNet.C(firstClu,secondClu) = costRecord(v).bridgeCost(i);
                superNet.C(secondClu,firstClu) = superNet.C(firstClu,secondClu);
            end
        end

        % endRoute update
        endClu = superRoute(v,superRouteL);
        superNet.C(end,endClu) = costRecord(v).residueCost(2);
        superNet.C(endClu,end) = superNet.C(end,endClu);
    end

    totalScoreHistory  = vertcat(totalScoreHistory, totalScore);
    totalScoreHistoryL = vertcat(totalScoreHistoryL, totalScoreL);
    totalTourHistory{trialNum} = tourRecord;
    superRouteHistory{trialNum} = superRoute;
    vehScoreHistory{trialNum} = vehScore;
    
    if ~firstTry %&& trialNum ~= 2
        if totalScoreHistory(end-1) < totalScore
            disp("Solution deteriorated.")
            disp("Termination condition met. Finishing the solver!");
            terminationType = "SOLDET";
            break;
        end
    end

%     if 1
%         terminationType = "HLPCONV";
%         break;
%     end


    firstTry = false;
    trialNum = trialNum + 1;
    prevSuperRoute = superRoute;
end
    
clusteringTime
interCompleteTime + intraCompleteTime
solveTime
totalScoreHistory
totalScoreHistoryL
vehScore

%% plot


[~,I] = min(totalScoreHistory);
finalTourRecord = totalTourHistory{I};

figure(4)
clf
grid on
hold on
for i = 1:size(G.Edges,1)
    startIdx = G.Edges.EndNodes(i,1);
    EndIdx = G.Edges.EndNodes(i,2);
    startPos = node(startIdx,:);
    EndPos = node(EndIdx,:);
    line([startPos(1) EndPos(1)],[startPos(2) EndPos(2)],[startPos(3) EndPos(3)]);
end
alpha = 0.1;
for i = 1:cluNum
    temp = node(nodeInCluIdx{i},:);
    plot3(temp(:,1),temp(:,2),temp(:,3),'x','LineWidth',5,'MarkerSize',5);
end
% mesh(voxelPosX,voxelPosY,voxelFilterData);
plot3(superNet.pos(:,1),superNet.pos(:,2),superNet.pos(:,3),'yx','MarkerSize',10,'LineWidth',5)
for i = 1:cluNum
    text(superNet.pos(i,1),superNet.pos(i,2),superNet.pos(i,3),num2str(i));
end
% for i = 1:N
%     text(node(i,1),node(i,2),0,num2str(i));
% end
hold on
for v= 1:vnum
    plot3(node(finalTourRecord{v},1),node(finalTourRecord{v},2),node(finalTourRecord{v},3)+2,'LineWidth',3)
end
axis equal

figure(6)
clf
hold on
grid on
plot(totalScoreHistory,'o-');
xlabel('trialNum');
ylabel('Score');
title('total score history')
%% post processing
% 2-opt
