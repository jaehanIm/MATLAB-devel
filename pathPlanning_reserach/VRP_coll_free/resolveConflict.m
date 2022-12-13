function [waitFee, reservation, reservationNum, occupancy, blocked, occu_hist, getAwayTime, unableFlag, damnFlag] = resolveConflict(reservation, reservationNum, occu_hist, vehTourLen, blocked, A, C, curNode, nextNode, tick, implicitRoute, vehIdx, getAwayTime)
unableFlag = false;
damnFlag = false;
global debugTemp

% get route info
if A(curNode, nextNode) == 0
    routeInfo = implicitRoute{curNode, nextNode};
else
    routeInfo = [curNode, nextNode];
end
L = size(routeInfo,2);
waitFee = 0;
occupancy = zeros(L-1,4);

% detect conflict & update
for i = 1:L-1 % for all route nodes
    initNode = routeInfo(i);
    termNode = routeInfo(i+1);
    reqInit = tick;
    reqTerm = tick + C(initNode,termNode);
    reqInfo = [reqInit, reqTerm];
    origTick = tick;
    
    % detect conflict and delay schedule
    reserveNum = reservationNum(initNode,termNode);
    runawayTime = min(getAwayTime{termNode}(getAwayTime{termNode}>origTick));
    runawayTimeInit = min(getAwayTime{initNode}(getAwayTime{initNode}>origTick));
    if isempty(runawayTime)
        runawayTime = inf;
    end
    if isempty(runawayTimeInit)
        runawayTimeInit = inf;
    end
    if reserveNum ~=0
        for j = 1:reserveNum % compare with for all reservations
            % conflict detection
            debugTemp.reserveNum_i = reserveNum;
            debugTemp.j_i = j;
            debugTemp.i_i = i;
            debugTemp.initNode_i = initNode;
            debugTemp.termNode_i = termNode;
            debugTemp.info_i = reservation{initNode,termNode}.info;
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
            if tick > runawayTimeInit && ~isempty(occu_hist) && i == 1
                damnFlag = true;
                break;
            end
        end
    end

    if damnFlag
        break;
    end

    if i<L-1 % check whether i can escape from the next node
        if ~isempty(runawayTime)
            preview_tick = reqTerm;
            preview_initNode = routeInfo(i+1);
            preview_termNode = routeInfo(i+2);
            preview_reqInit = preview_tick;
            preview_reqTerm = preview_reqInit + C(preview_initNode,preview_termNode);
            preview_reqInfo = [preview_reqInit,preview_reqTerm];
            preview_reserveNum = reservationNum(preview_initNode,preview_termNode);
            if preview_reserveNum ~= 0
                for pj = 1:preview_reserveNum
                    preview_reservationInfo = reservation{preview_initNode,preview_termNode}.info{pj};
                    preview_isConflict = detectConflict(preview_reqInfo,preview_reservationInfo);
                    if preview_isConflict
                        if preview_tick <= preview_reservationInfo(2)
                            preview_tick = preview_reservationInfo(2);
                        end
                    end
                    if preview_tick > runawayTime
                        unableFlag = true;
                        break;
                    end
                end
            end
        end
    end

    if unableFlag
        break;
    end

    if unableFlag == false
        % update waitFee for i (certain route node)
        localDelay = tick - origTick;
        waitFee = waitFee + localDelay;
        tick = reqTerm;
    
        occupancy(i,:) = [initNode, termNode, reqInit, reqTerm];
        getAwayTime{termNode} = vertcat(getAwayTime{termNode},reqInit);
    
        % if delay occurred -> extend previous reservation
        if localDelay > 0
            relevantNodes = find(A(initNode,:));
            for j = relevantNodes % extend initNode related reservations
                if ~isempty(reservation{j,initNode}.info)
                    localReserveLen = reservationNum(j,initNode);
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
                    localReserveLen = reservationNum(initNode,prevNode);
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
    %                     if reservation{initNode,prevSessionInit}.info{end-k+1}(3) == vehIdx
                        if reservation{initNode,prevSessionInit}.info{end-k+1}(2) == prevSessionEndingT && reservation{initNode,prevSessionInit}.info{end-k+1}(3) == vehIdx
                            updateIdx = localReserveLen-k+1;
                            break;
                        end
                    end
                    if ~isempty(updateIdx)
                        reservation{initNode,prevSessionInit}.info{updateIdx}(2) = reservation{initNode,prevSessionInit}.info{updateIdx}(2) + localDelay;
                        occu_hist{vehTourLen}(end,4) = reservation{initNode,prevSessionInit}.info{updateIdx}(2);
                    else
    %                     reservation{initNode,prevSessionInit}.info{end}(2) = reservation{initNode,prevSessionInit}.info{end}(2) + localDelay;
    %                     occu_hist{vehTourLen}(end,2) = reservation{initNode,prevSessionInit}.info{end}(2);
                    end
                end
            end
        end
        % update reservation schedule
        relevantNodes = find(A(termNode,:));
        for j = relevantNodes
            reservationNum(j,termNode) = reservationNum(j,termNode) + 1;
            reservation{j,termNode}.info{end+1} = [reqInit, reqTerm, vehIdx, localDelay];
        end
        reservationNum(termNode,initNode) = reservationNum(termNode,initNode) + 1;
        reservation{termNode,initNode}.info{end+1} = [reqInit, reqTerm, vehIdx, localDelay];
        if i == L-1
            blocked(nextNode) = 1;
            blocked(curNode) = 0;
        end
    
        % sort reservation
        for j = relevantNodes
            tempCell = [];
            tempCell{1,reservationNum(j,termNode)} = [];
            forSortInits = zeros(reservationNum(j,termNode),1);
            for k = 1:reservationNum(j,termNode)
                forSortInits(k) = reservation{j,termNode}.info{k}(1);
            end
            [~,I] = sort(forSortInits);
            for k = 1:reservationNum(j,termNode)
                tempCell{k} = reservation{j,termNode}.info{I==k};
            end
            reservation{j,termNode}.info = tempCell;
        end
        tempCell = [];
        tempCell{1,reservationNum(termNode,initNode)} = [];
        forSortInits = zeros(reservationNum(termNode,initNode),1);
        for k = 1:reservationNum(termNode,initNode)
            forSortInits(k) = reservation{termNode,initNode}.info{k}(1);
        end
        [~,I] = sort(forSortInits);
        for k = 1:reservationNum(termNode,initNode)
            tempCell{k} = reservation{termNode,initNode}.info{I==k};
        end
        reservation{termNode,initNode}.info = tempCell;
    else
        break;
    end
end

end