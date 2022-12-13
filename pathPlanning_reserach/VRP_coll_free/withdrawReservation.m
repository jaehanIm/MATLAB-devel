function [reservation, reservationNum] = withdrawReservation(occupancy, reservation, reservationNum, vehNum)
global debugTemp
occupancyLen = size(occupancy,1);

for i = 1:occupancyLen
    initNode = occupancy(i,1);
    termNode = occupancy(i,2);
    initTime = occupancy(i,3);
    index = [];
    debugTemp.i_w = i;
    debugTemp.occupancy_w = occupancy;
    debugTemp.vehNum_w = vehNum;
    debugTemp.initNode_w = initNode;
    debugTemp.termNode_w = termNode;
    for j = 1 : reservationNum(initNode,termNode)
        debugTemp.j_w = j;
        debugTemp.reservation_w = reservation{initNode,termNode}.info;
        if reservation{initNode,termNode}.info{j}(1) == initTime && reservation{initNode,termNode}.info{j}(3) == vehNum
            index = j;
            break;
        end
    end
    reservation{initNode,termNode}.info(index) = [];
    reservationNum(initNode,termNode) = reservationNum(initNode,termNode) -1;
end

end