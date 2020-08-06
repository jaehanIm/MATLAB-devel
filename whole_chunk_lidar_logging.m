% clear all

%% Variable setting

S1 = 7; %target focus distance [m]
N  = 3.5; % f
f  = 43; % focal length [mm]
S2 = 5:0.01:12;
S2 = S2 * 1000;
S1 = S1 * 1000;
coc = abs(S2-S1)./S2.*f^2./(N.*(S1-f));

%% Data loading
% logPlot
% flag = find(fix(data.jobSeq/max(data.jobSeq)),1);
% flage = flag + 1000;

% data = data(flag-100:flage,:); time = time(flag-100:flage);
yaw = data.rpy_2*pi/180;
lidar = [data.LidarRaw_0,data.LidarRaw_1,data.LidarRaw_2,data.LidarRaw_3,data.LidarRaw_4,data.LidarRaw_5,data.LidarRaw_6,data.LidarRaw_7];

pos = [data.posNed_0,data.posNed_1,data.posNed_2];
lidar_angle = [-21 -15 -9 -3 3 9 15 21] * pi/180;
lidar_primary_angle = data.bdLidarAngle * pi/180;
lidar_primary_dist = data.bdLidarDist;

for i = 1:size(lidar,1)
    for j = 1:size(lidar,2)
        if lidar(i,j) == -1
            lidar(i,j) = nan;
        end
    end
    if lidar_primary_dist(i) == -1
        lidar_primary_dist(i) = nan;
    end
end

d_x = []; d_y = [];
dp_x = []; dp_y = [];
for i = 1:length(lidar)
    d_x = vertcat(d_x,pos(i,1)+lidar(i,:).*cos(lidar_angle+yaw(i)));
    d_y = vertcat(d_y,pos(i,2)+lidar(i,:).*sin(lidar_angle+yaw(i)));
%     dp_x = vertcat(dp_x,pos(i,1)+lidar_primary_dist(i).*cos(lidar_primary_angle(i)+yaw(i)));
    dp_x = vertcat(dp_x,pos(i,1)+lidar_primary_dist(i).*cos(lidar_primary_angle(i)+yaw(i)));
    dp_y = vertcat(dp_y,pos(i,2)+lidar_primary_dist(i).*sin(lidar_primary_angle(i)+yaw(i)));
end
% 

%%Visualization
pos(:,3) = -pos(:,3);

% figure(1)
% clf
% hold on
% grid on
% % xlim([min(pos(:,1)),max(pos(:,1))]);
% % ylim([min(pos(:,2)),max(pos(:,2))]);
% plot3(pos(:,1),pos(:,2),pos(:,3),'k')
% plot3(pos(1,1),pos(1,2),pos(:,3),'r*','MarkerSize',5)
% plot3(pos(end,1),pos(end,2),pos(end,3),'bo','MarkerSize',5)
% 
% plot3(d_x,d_y,pos(:,3),'b.','MarkerSize',8)
% plot3(dp_x,dp_y,pos(:,3),'r.:','MarkerSize',5)
% xlabel('East [m]','FontSize',14)
% ylabel('North [m]','FontSize',14)
% 
% figure(2)
% clf
% hold on
% grid on
% xlim([-50 10]);
% ylim([-35 25]);
% an = animatedline('Marker','o');
% bd = animatedline('Marker','x');
% pn = animatedline('Marker','x','Color','r');
% plot(pos(:,1),pos(:,2),'k')
% plot(pos(1,1),pos(1,2),'r*','MarkerSize',5)
% plot(pos(end,1),pos(end,2),'bo','MarkerSize',5)
% line = []; line_2 = []; txt1 = []; txt2 = []; txt3 = [];
% line_3 = [];
% for i = 1:size(pos,1)
%     delete(line);
%     delete(line_2);
%     delete(line_3);
%     delete(txt1);delete(txt2); delete(txt3);
%     clearpoints(an)
%     clearpoints(bd)
%     clearpoints(pn)
%     addpoints(bd,pos(i,1),pos(i,2));
%     addpoints(an,d_x(i,1),d_y(i,1));
%     addpoints(an,d_x(i,2),d_y(i,2));
%     addpoints(an,d_x(i,3),d_y(i,3));
%     addpoints(an,d_x(i,4),d_y(i,4));
%     addpoints(an,d_x(i,5),d_y(i,5));
%     addpoints(an,d_x(i,6),d_y(i,6));
%     addpoints(an,d_x(i,7),d_y(i,7));
%     addpoints(an,d_x(i,8),d_y(i,8));
%     addpoints(pn,dp_x(i),dp_y(i));
%     line = plot_arrow([pos(i,1),pos(i,2),0],[pos(i,1)+10*cos(yaw(i)),pos(i,2)+10*sin(yaw(i)),0]);
%     [~,~,line_2] = draw_circle([pos(i,1),pos(i,2)],7,'k:');
%     line_3 = plot([pos(i,1) dp_x(i)],[pos(i,2) dp_y(i)],'r:');
%     txt1 = text(0,15,["time: " num2str(seconds(time(i)))]);
%     txt2 = text(-10,15,["jobseq: " num2str(data.jobSeq(i))]);
%     txt3 = text(-20,15,["lidar dist: " num2str(data.LidarDist(i))]);
%     drawnow;
%     pause(0.1)
% end

figure(3)
clf
hold on
grid on
xlim([min(pos(:,1))-50, max(pos(:,1))+50]);
ylim([min(pos(:,2))-50, max(pos(:,2))+50]);
zlim([20 80]);
an = animatedline('Marker','.');
bd = animatedline('Marker','x');
pn = animatedline('Marker','x','Color','r');
plot3(pos(:,1),pos(:,2),pos(:,3),'k')
plot3(pos(1,1),pos(1,2),pos(1,3),'r*','MarkerSize',5)
plot3(pos(end,1),pos(end,2),pos(end,3),'bo','MarkerSize',5)
% view(45,45);
line = []; line_2 = []; txt1 = []; txt2 = []; txt3 = [];
line_3 = [];
for i = 1:size(pos,1)
    delete(line);
    delete(line_2);
    delete(line_3);
    delete(txt1);delete(txt2); delete(txt3);
%     clearpoints(an) % lidar points
    clearpoints(bd)
    clearpoints(pn)
    addpoints(bd,pos(i,1),pos(i,2),pos(i,3));
    addpoints(an,d_x(i,1),d_y(i,1),pos(i,3));
    addpoints(an,d_x(i,2),d_y(i,2),pos(i,3));
    addpoints(an,d_x(i,3),d_y(i,3),pos(i,3));
    addpoints(an,d_x(i,4),d_y(i,4),pos(i,3));
    addpoints(an,d_x(i,5),d_y(i,5),pos(i,3));
    addpoints(an,d_x(i,6),d_y(i,6),pos(i,3));
    addpoints(an,d_x(i,7),d_y(i,7),pos(i,3));
    addpoints(an,d_x(i,8),d_y(i,8),pos(i,3));
    addpoints(pn,dp_x(i),dp_y(i),pos(i,3));
    line = plot_arrow([pos(i,1),pos(i,2),pos(i,3)],[pos(i,1)+10*cos(yaw(i)),pos(i,2)+10*sin(yaw(i)),pos(i,3)]);
    [~,~,line_2] = draw_circle3([pos(i,1),pos(i,2),pos(i,3)],7,'k:');
    line_3 = plot3([pos(i,1) dp_x(i)],[pos(i,2) dp_y(i)],[pos(i,3)  pos(i,3)],'r:');
    txt1 = text(0,15,["time: " num2str(seconds(time(i)))]);
    txt2 = text(-10,15,["jobseq: " num2str(data.jobSeq(i))]);
    txt3 = text(-20,15,["lidar dist: " num2str(data.LidarDist(i))]);
    drawnow;
%     pause(0.1)
end
