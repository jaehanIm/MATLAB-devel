function [waitFee, reservation, occupancy, blocked, prevDelay] = resolveConflict(reservation, occu_hist, vehTourLen, blocked, A, C, curNode, nextNode, tick, implicitRoute, vehIdx)
% get route info
if A(curNode, nextNode) == 0
    routeInfo = implicitRoute{curNode, nextNode};
else
    routeInfo = [curNode, nextNode];
end
L = size(routeInfo,2);
waitFee = 0;
occupancy = zeros(L-1,4);
prevDelay = 0;

% detect conflict & update
for i = 1:L-1 % for all route nodes
    initNode = routeInfo(i);
    termNode = routeInfo(i+1);
    reqInit = tick;
    reqTerm = tick + C(initNode,termNode);
    reqInfo = [reqInit, reqTerm];
    origTick = tick;
    localDelay = 0;
    
    % detect conflict and delay schedule
    reserveNum = reservation{initNode,termNode}.num;
    if reserveNum ~=0
        for j = 1:reserveNum % compare with for all reservations
            % conflict detection
            reservationInfo = reservation{initNode,termNode}.info{j};
            isConflict = detectConflict(reqInfo, reservationInfo);
            if isConflict % update reqinfo when conflict occurs
                if tick <= reservationInfo(2)
                    tick = reservationInfo(2);
                    reqInit = tick;
                    reqTerm = tick + C(initNode,termNode);
                    reqInfo = [reqInit,reqTerm];
                end
            end
        end
    end
    % update waitFee for i (certain route node)
    localDelay = tick - origTick;
    waitFee = waitFee + localDelay;
    tick = reqTerm;

    occupancy(i,:) = [initNode, termNode, reqInit, reqTerm];

%     if ~isempty(occu_hist)
%         occu_hist(1:vehTourLen);
%     end
%     colony.queen.occupancy{5,colony.queen.vehTourLen(5)}(end,:)

    % if delay occurred -> extend previous reservation
    if localDelay > 0
        relevantNodes = find(A(initNode,:));
        for j = relevantNodes % extend initNode related reservations
            if ~isempty(reservation{j,initNode}.info)
                localReserveLen = reservation{j,initNode}.num;
                updateIdx = [];
                for k = 1:localReserveLen
                    if reservation{j,initNode}.info{end-k+1}(3) == vehIdx %뒤부터 확인
                        updateIdx = localReserveLen-k+1;
                        break;
                    end
                end
                if ~isempty(updateIdx)
                    reservation{j,initNode}.info{updateIdx}(2) = reservation{j,initNode}.info{updateIdx}(2) + localDelay;
                else
%                     reservation{j,initNode}.info{end}(2) = reservation{j,initNode}.info{end}(2) + localDelay;
                end
%                 
            end
        end
        % extend from init to prev reservation
        if i > 1 % if we are within current session
            prevNode = routeInfo(i-1);
            if ~isempty(reservation{initNode,prevNode}.info)
                localReserveLen = reservation{initNode,prevNode}.num;
                updateIdx = [];
                for k = 1:localReserveLen
                    if reservation{initNode,prevNode}.info{end-k+1}(3) == vehIdx %뒤부터 확인
                        updateIdx = localReserveLen-k+1;
                        break;
                    end
                end
                if ~isempty(updateIdx)
                    reservation{initNode,prevNode}.info{updateIdx}(2) = reservation{initNode,prevNode}.info{updateIdx}(2) + localDelay;
                else
%                     reservation{initNode,prevNode}.info{end}(2) = reservation{initNode,prevNode}.info{end}(2) + localDelay;
                end
            end
        else % if we are extending previous session
            if ~isempty(occu_hist)
                prevSessionInit = occu_hist{vehTourLen}(end,1);
                prevSessionTerm = occu_hist{vehTourLen}(end,2);
                prevSessionEndingT = occu_hist{vehTourLen}(end,4);
                localReserveLen = size(reservation{initNode,prevSessionInit}.info,2);
                updateIdx = [];
                for k = 1:localReserveLen
                    if reservation{initNode,prevSessionInit}.info{end-k+1}(3) == vehIdx
%                     if reservation{initNode,prevSessionInit}.info{end-k+1}(2) == prevSessionEndingT && reservation{initNode,prevSessionInit}.info{end-k+1}(3) == vehIdx
                        updateIdx = localReserveLen-k+1;
                        break;
                    end
                end
                if ~isempty(updateIdx)
                    reservation{initNode,prevSessionInit}.info{updateIdx}(2) = reservation{initNode,prevSessionInit}.info{updateIdx}(2) + localDelay;
                else
%                     reservation{initNode,prevSessionInit}.info{end}(2) = reservation{initNode,prevSessionInit}.info{end}(2) + localDelay;
                end
            end
        end
    end

    % update reservation schedule
    relevantNodes = find(A(termNode,:));
    for j = relevantNodes
        reservation{j,termNode}.num = reservation{j,termNode}.num + 1;
        reservation{j,termNode}.info{end+1} = [reqInit, reqTerm, vehIdx, localDelay];
    end
    reservation{termNode,initNode}.num = reservation{termNode,initNode}.num + 1;
    reservation{termNode,initNode}.info{end+1} = [reqInit, reqTerm, vehIdx, localDelay];
    if i == L-1
        blocked(nextNode) = 1;
        blocked(curNode) = 0;
    end

    % sort reservation
    for j = relevantNodes
        tempCell = [];
        tempCell{1,reservation{j,termNode}.num} = [];
        forSortInits = zeros(reservation{j,termNode}.num,1);
        for k = 1:reservation{j,termNode}.num
            forSortInits(k) = reservation{j,termNode}.info{k}(1);
        end
        [~,I] = sort(forSortInits);
        for k = 1:reservation{j,termNode}.num
            tempCell{k} = reservation{j,termNode}.info{I==k};
        end
        reservation{j,termNode}.info = tempCell;
    end
    tempCell = [];
    tempCell{1,reservation{termNode,initNode}.num} = [];
    forSortInits = zeros(reservation{termNode,initNode}.num,1);
    for k = 1:reservation{termNode,initNode}.num
        forSortInits(k) = reservation{termNode,initNode}.info{k}(1);
    end
    [~,I] = sort(forSortInits);
    for k = 1:reservation{termNode,initNode}.num
        tempCell{k} = reservation{termNode,initNode}.info{I==k};
    end
    reservation{termNode,initNode}.info = tempCell;
end

end