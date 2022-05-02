function out = leastTourAgent(tourInfo,L,vehNum,vehTourLen)
cumLen = zeros(vehNum,1);
for i = 1:vehNum
    tour = tourInfo;
    for j = 1:vehTourLen(i)-1
        cumLen(i,1) = cumLen(i,1) + L(tour(i,j),tour(i,j+1));
    end
end
[~,out]=min(cumLen);
end