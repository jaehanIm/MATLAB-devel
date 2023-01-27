tour = colony.queen.tour;
totalTimeHistory{vnum} = [];
totTourLen = colony.queen.vehTourLen;
detailRoute = {};
detailRoute{vnum} = [];

%% Decipher tour
for i = 1:vnum
    tourLen = length(tour(i,tour(i,:)~=0));
    for j = 1:tourLen-2
        start = tour(i,j);
        finish = tour(i,j+1);
        detailRoute{i} = vertcat(detailRoute{i},start);
        if A(start,finish) == 0
            interLen = length(implicitRoute{start,finish})-2;
            for k = 1:interLen
                detailRoute{i} = vertcat(detailRoute{i},implicitRoute{start,finish}(k+1));
            end
        end
    end
    detailRoute{i} = vertcat(detailRoute{i},1);
end

%% asdf
tour = detailRoute;
for i = 1:vnum
    tourLen = length(tour{i});
    timeHistory = zeros(1,tourLen);
    timeHistory(1) = 0;
    for j = 2:tourLen
        timeHistory(j) = timeHistory(j-1) + C(tour{i}(j-1),tour{i}(j));
    end
    totalTimeHistory{i}  = timeHistory;
end

soleOccup{N,N}.num = 0;
soleOccup{N,N}.info = [];
for i = 1:N
    for j = 1:N
        soleOccup{i,j}.num = 0;
        soleOccup{i,j}.info = [];
    end
end

for i = 1:vnum
    tourLen = length(tour{i});
    for j = 1:tourLen-1
        start = tour{i}(j);
        finish = tour{i}(j+1);
        startT = totalTimeHistory{i}(j);
        finishT = totalTimeHistory{i}(j+1);
        soleOccup{start,finish}.num = soleOccup{start,finish}.num+1;
        soleOccup{finish,start}.num = soleOccup{finish,start}.num+1;
        soleOccup{start,finish}.info = vertcat(soleOccup{start,finish}.info,[startT, finishT, i]);
        soleOccup{finish, start}.info = vertcat(soleOccup{finish,start}.info,[startT, finishT, i]);
    end
end

warn = 0;
accumHist = [];
for i = 1:N-1
    for j = i+1:N
        localSchedule = soleOccup{i,j}.info;
        accum = checkConflict(localSchedule);
        accumHist = vertcat(accumHist,accum);x
        warn = warn + accum;
    end
end

function accum = checkConflict(info)
accum = 0;
for i = 1:size(info,1)-1
    for j = i+1:size(info,1)
        first = info(i,1:2);
        firstVeh = info(i,3);
        second = info(j,1:2);
        secondVeh = info(j,3);
        if detectConflict(first,second) && firstVeh ~= secondVeh
            accum = accum + min([first(2),second(2)]) - max([first(1),second(1)]);
        end
    end
end

end
