data1 = load('oksang_1.mat');
data1_vel = data1.meanVelBar;
data1_std = data1.stdVelBar;
data1_tot = data1.meanTotVel;
data1_totStd = data1.stdTotVel;
data2 = load('oksang_2.mat');
data2_vel = data2.meanVelBar;
data2_std = data2.stdVelBar;
data2_tot = data2.meanTotVel;
data2_totStd = data2.stdTotVel;
data4 = load('oksang_3.mat');
data4_vel = data4.meanVelBar;
data4_std = data4.stdVelBar;
data4_tot = data4.meanTotVel;
data4_totStd = data4.stdTotVel;

totVelData = [data1_vel,data2_vel,data4_vel];
totStdData = [data1_std;data2_std;data4_std]';
totTotVelData = [data1_tot;data2_tot;data4_tot];
totStdTotVelData = [data1_totStd;data2_totStd;data4_totStd];

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
errorbar(x',totVelData,totStdData,'k.','LineWidth',1)
grid on
ylabel('Pointing Velocity [m/s]')
legend(b,{'Original','Updated','Vehicle'})
ylim([0 inf])
title('Average pointing vel.')

figure(2)
clf
bar(1,totTotVelData(1))
% xticks([1])
hold on
bar(2,totTotVelData(2))
bar(3,totTotVelData(3))
xticks([1,2,3]);
xticklabels({'Original','Updated','Vehicle'});
errorbar(totTotVelData,totStdTotVelData,'k.','LineWidth',1)
ylabel('Pointing Velocity [m/s]')
title('Average pointing vel. (Hovering)')
grid on