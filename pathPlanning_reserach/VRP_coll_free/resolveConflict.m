function [waitFee, reservation, occupancy] = resolveConflict(reservation, A, C, curNode, nextNode, tick, implicitRoute, vehIdx)
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
    localDelay = 0;
    
    % detect conflict and delay schedule
    reserveNum = reservation{initNode,termNode}.num;
    if reserveNum ~=0
        for j = 1:reserveNum % compare with for all reservations
            % conflict detection
            reservationInfo = reservation{initNode,termNode}.info{j};
            isConflict = detectConflict(reqInfo, reservationInfo);
            if isConflict % update reqinfo when conflict occurs
                if tick < reservationInfo(2)
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

    % update reservation schedule
    relevantNodes = find(A(termNode,:));
    for j = relevantNodes
%         reservation{termNode,j}.num = reservation{termNode,j}.num + 1;
        reservation{j,termNode}.num = reservation{j,termNode}.num + 1;
%         reservation{termNode,j}.info{end+1} = [reqInit, reqTerm, vehIdx, localDelay];
        reservation{j,termNode}.info{end+1} = [reqInit, reqTerm, vehIdx, localDelay];
    end
end

% sort reservation schedule (according to start time)


end