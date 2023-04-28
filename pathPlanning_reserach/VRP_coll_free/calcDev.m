function devSum = calcDev(tour, implicitRoute, curNode, nextNode, nodePos, C, egoVIdx, N)
% tour : NxtourLen
% implicitRoute : NxN
% nodePos : Nx3
% C : NxN cost matrix
% output : 1xN matrix

global resolution

vnum = size(tour,1);
devSum = 0;
startTime = getTimeAtArrival(tour(egoVIdx,:), C);
arrivTime = startTime + C(curNode,nextNode);
candidateRoute = nonzeros(tour(egoVIdx,:));
candidateRoute = [candidateRoute;nextNode];

detailTour = cell(vnum,1);
for v = 1:vnum
    detailTour{v} = routeDecoder(nonzeros(tour(v,:)),implicitRoute, N);
end
detailCandidateRoute = routeDecoder(nonzeros(candidateRoute),implicitRoute,N);



for v = 1:vnum
    if v ~= egoVIdx
        for time = linspace(startTime , arrivTime, resolution)
            egoPos = getVehPos(detailCandidateRoute, time, C, nodePos, implicitRoute, N);
            oppPos = getVehPos(detailTour{v}, time, C, nodePos, implicitRoute, N);
            dev = norm(egoPos - oppPos);
            devSum = devSum + dev;
        end
    end
end

end