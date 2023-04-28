function pos = getVehPos(totTour, time, C, nodePos, implicitRoute, N)
% nonzeros filtering required!!! before using this function
% tour = routeDecoder(totTour, implicitRoute, N);
% tour = nonzeros(totTour);

tour = totTour;
tourLen = length(tour)-1;
curTfix = 0;
nextTfix = 0;
isUnderTour = false;

if tourLen ~= 0
    for i = 1:tourLen
        curTfix = nextTfix;
        curNode = tour(i);
        nextNode = tour(i+1);
        nextTfix = curTfix + C(curNode,nextNode);
        if nextTfix > time
            isUnderTour = true;
            break;
            % check whether this part is read
        end
    end
    if isUnderTour == true
        a = time - curTfix;
        c = nextTfix - curTfix;
        curPos = nodePos(curNode,:);
        nextPos = nodePos(nextNode,:);
        pos = curPos + (nextPos - curPos) * a / c;
    else
        pos = nodePos(tour(end-1),:);
    end   
else
    % No tour yet
    pos = nodePos(1,:);
end

end