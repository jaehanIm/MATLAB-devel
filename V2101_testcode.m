[data,time] = loader("/home/jaehan/log/gdLogCsv/gdLog_201216_194606.csv");
A = [[-3.25884 46.1553 -33.1394]; [-1.49107 47.923 -37.4695]; [2.04446 44.3875 -37.4695]; [0.276692 42.6197 -33.1394]];
A = vertcat(A,A(1,:));

doingJobIdx = find(data.fcMcMode == 2);
startTime = doingJobIdx(find(doingJobIdx>600,1));

l = 25;
d2r = pi/180;
r2d = 1/d2r;

gimbalRPYCmd_0 = data.gimbalRPYCmd_0(startTime:end)*d2r;
gimbalRPYCmd_1 = data.gimbalRPYCmd_1(startTime:end)*d2r;
gimbalRPYCmd_2 = data.gimbalRPYCmd_2(startTime:end)*d2r;
gimbalRPY_0 = data.gimbalRPY_0(startTime:end);
gimbalRPY_1 = data.gimbalRPY_1(startTime:end);
gimbalRPY_2 = data.gimbalRPY_2(startTime:end);
posNed_0 = data.posNed_0(startTime:end);
posNed_1 = data.posNed_1(startTime:end);
posNed_2 = data.posNed_2(startTime:end);
eulerY = data.rpy_2(startTime:end)*d2r;

dataSize = length(gimbalRPY_0);
WTdir = 30;
inspectionAngle = -30;
inspectionAngle = inspectionAngle * d2r;
WTdir = WTdir * d2r;

dcmI2ImageFrame = angle2dcm(WTdir,inspectionAngle,0,'zyx');
rotmI2gridFrame = angle2dcm(-gimbalRPYCmd_1,-gimbalRPYCmd_0,-gimbalRPYCmd_2-eulerY,'yxz');
rotmI2ImageFrame = angle2dcm(0,-inspectionAngle,0,'xyz');
% rotmI2ImageFrame = angle2dcm(-inspectionAngle,-30*d2r,-WTdir,'yxz');

targetImageFrame = [l,1.759/2,1.296/2;l,1.759/2,-1.296/2;l,-1.759/2,-1.296/2;l,-1.759/2,1.296/2;l,1.759/2,1.296/2];
targetInspectionFrame = [l,40,20;l,40,-20;l,-40,-20;l,-40,20;l,40,20];
targetNed = dcmI2ImageFrame' * targetImageFrame';
targetInsNed = rotmI2ImageFrame * targetInspectionFrame';
targetNed = targetNed' + [posNed_0(1),posNed_1(1),posNed_2(1)];
targetInsNed = targetInsNed' + [posNed_0(1),posNed_1(1),posNed_2(1)];

% draw sample
pointingNed = zeros(3,length(gimbalRPY_0));
for i = 1:length(gimbalRPY_0)
    pointingNed(:,i) = rotmI2gridFrame(:,:,i) * [l,0,0]';
end
pointingNed = pointingNed' + [posNed_0(1),posNed_1(1),posNed_2(1)];

% gridRectangles
gridRectangles = zeros(5,3,dataSize);
for i = 1:dataSize
    temp = rotmI2gridFrame(:,:,i) * targetImageFrame';
    
    gridRectangles(:,:,i) = temp' + [posNed_0(1),posNed_1(1),posNed_2(1)];
end

figure(1)
clf
hold on
grid on
plot(gimbalRPYCmd_2/d2r,gimbalRPYCmd_1/d2r)
axis equal
xlabel('Yaw angle')
ylabel('Pitch angle')

figure(2)
clf
subplot(2,2,2)
title('lateral view')
hold on
grid on
plot3(posNed_1(1),posNed_0(1),-posNed_2(1),'x','MarkerSize',10)
plot3(targetNed(:,2),targetNed(:,1),-targetNed(:,3),'k')
plot3(pointingNed(:,2),pointingNed(:,1),-pointingNed(:,3))
plot_arrow([posNed_1(1) posNed_0(1) -posNed_2(1)],[posNed_1(1) posNed_0(1) -posNed_2(1)] + l*[sin(eulerY(end)) cos(eulerY(end)) 0]/1.2);
view(45,45)
axis equal
view(90,0)

subplot(2,2,1)
hold on
grid on
title('3d view')
plot3(posNed_1(1),posNed_0(1),-posNed_2(1),'x','MarkerSize',10)
plot3(targetNed(:,2),targetNed(:,1),-targetNed(:,3),'k--')
plot3(pointingNed(:,2),pointingNed(:,1),-pointingNed(:,3))
plot_arrow([posNed_1(1) posNed_0(1) -posNed_2(1)],[posNed_1(1) posNed_0(1) -posNed_2(1)] + l*[sin(eulerY(end)) cos(eulerY(end)) 0]/1.2);
plot3(targetInsNed(:,2),targetInsNed(:,1),-targetInsNed(:,3),'r')
for i = 1:dataSize
    plot3(gridRectangles(:,2,i),gridRectangles(:,1,i),-gridRectangles(:,3,i),'k')
end
view(45,45)
axis equal

subplot(2,2,3)
hold on
grid on
title('birds eye view')
plot3(posNed_1(1),posNed_0(1),-posNed_2(1),'x','MarkerSize',10)
plot3(targetNed(:,2),targetNed(:,1),-targetNed(:,3),'k')
plot3(pointingNed(:,2),pointingNed(:,1),-pointingNed(:,3))
plot_arrow([posNed_1(1) posNed_0(1) -posNed_2(1)],[posNed_1(1) posNed_0(1) -posNed_2(1)] + l*[sin(eulerY(end)) cos(eulerY(end)) 0]/1.2);
view(45,45)
axis equal
view(0,90)

subplot(2,2,4)
hold on
grid on
plot(gimbalRPYCmd_0*r2d)
plot(gimbalRPYCmd_1*r2d)
plot(gimbalRPYCmd_2*r2d)
ylabel('Angle [deg]')
title('Gimbal Angle')
legend('Roll','Pitch','Yaw')
% 
figure(3)
clf
hold on
grid on
title('3d view')
plot3(posNed_1(1),posNed_0(1),-posNed_2(1),'x','MarkerSize',10)
plot3(targetNed(:,2),targetNed(:,1),-targetNed(:,3),'k--')
plot3(pointingNed(:,2),pointingNed(:,1),-pointingNed(:,3))
plot_arrow([posNed_1(1) posNed_0(1) -posNed_2(1)],[posNed_1(1) posNed_0(1) -posNed_2(1)] + l*[sin(eulerY(end)) cos(eulerY(end)) 0]/1.2);
plot3(targetInsNed(:,2),targetInsNed(:,1),-targetInsNed(:,3),'r')
for i = 1:dataSize
    plot3(gridRectangles(:,2,i),gridRectangles(:,1,i),-gridRectangles(:,3,i),'k')
end
view(45,45)
axis equal