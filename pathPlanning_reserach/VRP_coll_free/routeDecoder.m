function detail = routeDecoder(upperTour, implicitRoute,N)
% upperTour Nx1
% implicitRoute NxN
tour = zeros(N,1);
upperTourLen = length(upperTour)-1;
idxCount = 0;
for i = 1:upperTourLen    
    initNode = upperTour(i);
    finNode = upperTour(i+1);
    if ~isempty(implicitRoute{initNode,finNode})
        implList = implicitRoute{initNode,finNode}; implList(end) = [];
        implLen = length(implList);
        tour(idxCount+1:idxCount+implLen) = implList;
        idxCount = idxCount + implLen;
    else
        tour(idxCount+1) = initNode;
        idxCount = idxCount + 1;
    end
end
tour(idxCount+1) = upperTour(end);
detail = nonzeros(tour);
end