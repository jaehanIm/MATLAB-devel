yaw = -40*pi/180;
dcmI2WtTail = angle2dcm(yaw, 0, 0,'zyx');
rot_data = dcmI2WtTail * [data.posNed_1';data.posNed_0';data.posNed_2'];
% range = 4600:8441;
% range = 8440:9678;
% range = 9678:13170;
range = 14190:19340;

% figure(3)
% t = seconds(time);
% clf
% hold on
% grid on
% xlim([-10 60])
% ylim([-10 60])
% dot = animatedline('Marker','.');
% txt = []; txt2 = [];
% for i = 1:length(rot_data)
% delete(txt);
% delete(txt2);
% addpoints(dot,rot_data(1,i),rot_data(2,i));
% txt = text(30,10,['Time: ',num2str(t(i)),'s']);
% txt2 = text(30,5,['Time idx: ',num2str(i)]);
% drawnow
% end

figure(5)
clf
plot(time(range),rot_data(1,range)-mean(rot_data(1,range)))
hold on
plot(time(range),data.LidarDist(range)-mean(data.LidarDist(range)))
legend('position_','lidardist')
grid on
title('거리유지 pos & lidar')
xlabel('time')
ylabel('m')
text(130,0.4,'mean lidar dist : 11.04m')


% plot(dist1.coeff.breaks,ppval(dist1.coeff,dist1.coeff.breaks))
