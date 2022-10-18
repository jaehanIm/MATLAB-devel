function [ colony, incompleteFlag ] = createColonyVRPCF( graph, antNo, tau, eta, param, vehNum, capacity, servTime, implicitRoute)

global homeIdx debugTemp tau0
incompleteFlag = false;

alpha = param.alpha;
beta = param.beta;
% gamma = param.gamma;
q = param.q;
% lambda = param.lambda;
psi = param.psi;

A = graph.A;
N = graph.n;
C = graph.edges;
debugTemp.N = N;
colony = [];

for i = 1 : antNo
    
    unvisitedNum = N;
    
    initial_node = homeIdx;
    colony.ant(i).tour(:,1) = ones(vehNum,1) * initial_node;
    colony.ant(i).tick = zeros(vehNum,1);
    colony.ant(i).tickHistory = zeros(vehNum,1);
    colony.ant(i).occupancy = {};
    unvisitedNum = unvisitedNum - 1;
    vehTourLen = ones(vehNum,1);
    reservation{N,N} = [];
    blocked = zeros(N,1);
    stuckVeh = zeros(vehNum,1);
    getAwayTime = {};
    tempBlocked = {};

    for m = 1:N
        for n = 1:N
            reservation{m,n}.num = 0;
            reservation{m,n}.info = [];
        end
    end
    getAwayTime{N} = [];
    tempBlocked{vehNum} = [];

    while unvisitedNum ~= 0
            j = getActiveVeh(colony.ant(i).tick, stuckVeh);

            if unvisitedNum == 0
                break;
            end            
            
            % i = antNo, j = vehNum
            currentNode = colony.ant(i).tour(j,vehTourLen(j));
            P_allNodes = tau(currentNode,:).^alpha.*eta(currentNode,:).^beta;
            P_allNodes(nonzeros(colony.ant(i).tour(:))) = 0;
            P_allNodes(tempBlocked{j}) = 0;

            possibleNodes = P_allNodes~=0;
            debugTemp.possibleNodes2 = possibleNodes;
            for n = 1:N
                if possibleNodes(n) == true
                    % check if the node is really viable (not in blockedList)
                    if A(currentNode,n) == 0
                        localRoute = implicitRoute{currentNode,n};
                        localLen = length(localRoute);
                        for nn = 2:localLen
                            if blocked(localRoute(nn)) == 1
                                possibleNodes(n) = false;
                                break;
                            end
                        end
                    else
                        if blocked(n) == 1
                            possibleNodes(n) = false;
                        end
                    end
                end
            end

            debugTemp.possibleNodes = possibleNodes;
            debugTemp.blocked = blocked;

            P_allNodes(~possibleNodes) = 0;

            if ~isempty(nonzeros(P_allNodes)) % normal state
                % reset stuckVeh
                stuckVeh = zeros(vehNum,1);

                % Calculate pseudo-random P
                P = P_allNodes./sum(P_allNodes);
                
                % for debugging
                debugTemp.P_allNodes = P_allNodes;
                debugTemp.P = P;
                debugTemp.currentNode = currentNode;
                debugTemp.i = i;
                debugTemp.j = j;
                debugTemp.colony = colony;
                debugTemp.vehTourLen = vehTourLen;
                debugTemp.unvisitedNum = unvisitedNum;
                debugTemp.tempBlocked =tempBlocked;
                debugTemp.reservation = reservation;
                if unvisitedNum == N-1
                    debugTemp.initColony = colony;
                end
                
                if rand(1) > q
                    nextNode = rouletteWheel(P);
                else
                    [~,nextNode] = max(P_allNodes);
                    nextNode = nextNode(1); % just in case multiple nextNodes are generated
                end

                if vehTourLen(j) == 1%isempty(colony.ant(i).occupancy)
                    occu_hist = [];
                else
                    occu_hist = colony.ant(i).occupancy(j,:);
                    orig_occu_hist = occu_hist;
                end

                [timeSlack, reservation_o, occupancy_o, blocked_o, occu_hist_o, getAwayTime_o, unableFlag, damnFlag] = resolveConflict(reservation, occu_hist, vehTourLen(j), blocked, A, C, currentNode, nextNode, colony.ant(i).tick(j), implicitRoute, j, getAwayTime);

                if damnFlag && vehTourLen(j) > 1
                    disp("[notice] Withdraw!")
                    reservation = withdrawReservation(colony.ant(i).occupancy{j,vehTourLen(j)},reservation,j);
                    colony.ant(i).tour(j,vehTourLen(j)) = 0;
                    unvisitedNum = unvisitedNum+1;
                    colony.ant(i).occupancy{j,vehTourLen(j)} = [];
                    colony.ant(i).tick(j) = colony.ant(i).tickHistory(j,vehTourLen(j)-1);
                    colony.ant(i).tickHistory(j,vehTourLen(j)) = 0;
                    vehTourLen(j) = vehTourLen(j)-1;
%                     stuckVeh(j) = true;
                else
                    if ~unableFlag && ~damnFlag
                        tempBlocked{j} = [];
                        %update info
                        reservation = reservation_o;
                        blocked = blocked_o;
                        occupancy = occupancy_o;
                        occu_hist = occu_hist_o;
                        getAwayTime = getAwayTime_o;
        
                        if ~isempty(occu_hist)
                            colony.ant(i).occupancy{j,vehTourLen(j)} = occu_hist{vehTourLen(j)};
                        end
        
                        vehTourLen(j) = vehTourLen(j) +1;
                        unvisitedNum = unvisitedNum - 1;
                        colony.ant(i).tour(j,vehTourLen(j)) = nextNode;
                        colony.ant(i).occupancy{j,vehTourLen(j)} = occupancy;
                        colony.ant(i).tick(j) = colony.ant(i).tick(j) + C(currentNode,nextNode) + timeSlack;
                        colony.ant(i).tickHistory(j,vehTourLen(j)) = colony.ant(i).tick(j);
            
                        % local pheromone update
                        tau(currentNode,nextNode) = tau(currentNode,nextNode) * (1-psi) + tau0 * psi;
                    else
                        tempBlocked{j} = vertcat(tempBlocked{j},nextNode);
                        tempBlocked{j} = unique(tempBlocked{j});
                    end
                end
            else % no viable nodes to visit
                stuckVeh(j) = true;
                debugTemp.stuck = stuckVeh;
                if stuckVeh == true
                    disp("Shit...")
                    incompleteFlag = true;
                    break;
                end
            end
    end
    
    colony.ant(i).vehTourLen = vehTourLen;
    colony.ant(i).reservation = reservation;
    % complete the tour 
    for j = 1:vehNum
        colony.ant(i).tour(j,vehTourLen(j)+1) = homeIdx;
        colony.ant(i).tickHistory(j,vehTourLen(j)+1) = ...
            colony.ant(i).tickHistory(j,vehTourLen(j)) + C(colony.ant(i).tour(j,vehTourLen(j)),homeIdx);
        colony.ant(i).occupancy{j,vehTourLen(j)+1} = [0 0 0 0];
    end
end

end