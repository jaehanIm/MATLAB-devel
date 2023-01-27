% dataList = ["oksang_1.mat","oksang_2.mat","oksang_4.mat","oksang_3.mat"];
dataList = ["oksang_wind_1","oksang_wind_2","oksang_wind_3"];
dataNum = length(dataList);
legendArray = {'Neutral','Mid wind','High wind'};

data = {};
for i = 1:dataNum
    data{i}.data = load(dataList(i));
    data{i}.vel = data{i}.data.meanVelBar;
    data{i}.std = data{i}.data.stdVelBar;
    data{i}.tot = data{i}.data.meanTotVel;
    data{i}.totStd = data{i}.data.stdTotVel;
end

totVelData = [];
totStdData = [];
totTotVelData = [];
totStdTotVelData = [];

for i = 1:dataNum
    totVelData = horzcat(totVelData, data{i}.vel);
    totStdData = vertcat(totStdData, data{i}.std);
    totTotVelData = horzcat(totTotVelData, data{i}.tot);
    totStdTotVelData = horzcat(totStdTotVelData, data{i}.totStd);
end

figure(1)
clf
b = bar(totVelData, 'grouped');
hold on
xticks([1,2,3,4,5,6,7]);
xticklabels({'1-5hz','5-30hz','30-75hz','75-150hz','150-800hz','800-1200hz','1200~hz'});
[ngroups,nbars] = size(totVelData);
x = nan(nbars,ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
errorbar(x,totVelData',totStdData,'k.','LineWidth',1)
grid on
ylabel('Pointing Velocity [m/s]')
% legend(b,{'C','C++','D','Vehicle','Vehicle2'})
legend(b,legendArray)
ylim([0 inf])
title('Average pointing vel.')

figure(2)
clf
for i = 1:dataNum
    bar(i,totTotVelData(i));
    hold on
end
xticks([1,2,3,4]);
xticklabels(legendArray);
errorbar(totTotVelData,totStdTotVelData,'k.','LineWidth',1)
ylabel('Pointing Velocity [m/s]')
title('Average pointing vel. (Hovering)')
grid on