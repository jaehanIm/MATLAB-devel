clear all

d2r = pi/180;

[data,time] = loader('/home/jaehan/Desktop/test flight/200715_154607/gdLog_200715_154607.csv');
posNed = [data.posNed_0,data.posNed_1,data.posNed_2];
lidarData = data.LidarDist;

dcmI2bridge = angle2dcm(wrapToPi(40.52*d2r), 0, 0,'zyx');
posData = dcmI2bridge * posNed';
posData = posData(2,:);
posData = posData';

sample = [4600 8440]; % 1st lidar
% sample = [8440 9678]; % 1st hold
% sample = [9678 13170]; % 2nd lidar
% sample = [14190 19340]; % 3rd zlidar

figure(1)
clf
plot(time,posData);
grid on
title('Position history')

figure(2)
clf
plot(time,lidarData);
grid on
title('Lidar distance history')

figure(3)
clf
hold on
plot(time(sample(1):sample(2)),posData(sample(1):sample(2))-mean(posData(sample(1):sample(2))));
plot(time(sample(1):sample(2)),lidarData(sample(1):sample(2))-mean(lidarData(sample(1):sample(2))));
grid on
title('Section history')
legend('drone position','lidar distance')