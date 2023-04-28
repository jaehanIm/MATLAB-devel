function time = getTimeAtArrival(tourRaw, C)
% individual tour
tour = nonzeros(tourRaw);
L = length(tour)-1;
time = 0;

if L ~= 0
    for i = 1:L
        curNode = tour(i);
        nextNode = tour(i+1);
        time = time + C(curNode, nextNode);
    end
end

end