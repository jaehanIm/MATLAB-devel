clear all
% close all
d2r = pi/180;
shootDistance = 25;


spdpath = '/home/jaehan/Desktop/test_flight/0420 FT/210420_121338/spdlog_210420_121338.txt';
gdpath = '/home/jaehan/Desktop/test_flight/0420 FT/210420_121338/gdLog_210420_121338.csv';


fileID = fopen(spdpath);
[ftdata,time] = loader(gdpath);

n = 1;
while ~feof(fileID)
    data{n} = fgetl(fileID);
    n = n+1;
end

gimbalrpy_0 = ftdata.gimbalRPYCmd_0*d2r;
gimbalrpy_1 = ftdata.gimbalRPYCmd_1*d2r;
gimbalrpy_2 = ftdata.gimbalRPYCmd_2*d2r;
rpyCmd_2 = ftdata.ySp*d2r;
posNed_0 = ftdata.posCmdNed_0;
posNed_1 = ftdata.posCmdNed_1;
posNed_2 = ftdata.posCmdNed_2;

% extract interest point
isInterest = gimbalrpy_1 == 0;
gimbalrpy_0(isInterest) = [];
gimbalrpy_1(isInterest) = [];
gimbalrpy_2(isInterest) = [];
rpyCmd_2(isInterest) = [];
posNed_0(isInterest) = [];
posNed_1(isInterest) = [];
posNed_2(isInterest) = [];
ftsize = size(posNed_0,1);

% extract unique points
filtered = 0;
for i = 2:ftsize
    idx = i - filtered;
    if gimbalrpy_1(idx-1) == gimbalrpy_1(idx) && gimbalrpy_2(idx-1) == gimbalrpy_2(idx)
        gimbalrpy_0(idx) = [];
        gimbalrpy_1(idx) = [];
        gimbalrpy_2(idx) = [];
        rpyCmd_2(idx) = [];
        posNed_0(idx) = [];
        posNed_1(idx) = [];
        posNed_2(idx) = [];        
        filtered = filtered + 1;
    end
end

% Extract wtYaw angle
wtYaw = [];
for i = 1:size(data,2)
is_wtYaw = strfind(data{i}, 'WT yaw angle and current vehicle yaw');
    if ~isempty(is_wtYaw)
        wtYaw = str2double(data{i}(is_wtYaw+38:is_wtYaw+44));
    end
end
wtYaw = -wtYaw + pi;


% Calculate DCM
ftsize = length(posNed_0);
posNed = [posNed_0,posNed_1,posNed_2]';

dcmI2ShootSurface = angle2dcm(wtYaw+pi/2,0*pi/180,0,'zyx'); % #0 LE
% dcmI2ShootSurface = angle2dcm(wtYaw-pi/2,-60*pi/180,0,'zyx'); % #0 LE
dcmI2ShootSurface_2 = dcmI2ShootSurface;
dcmI2Gimbal = angle2dcm(gimbalrpy_2,gimbalrpy_0,gimbalrpy_1,'zxy');
dcmI2Nav = angle2dcm(zeros(ftsize,1),zeros(ftsize,1),rpyCmd_2,'xyz');

maxvertexN = [];
maxvertexE = [];
maxvertexD = [];
minvertexN = [];
minvertexE = [];
minvertexD = [];
summary = [];

gridPosImageFrame = [];
gridPosNED = [];
addedVertexImageFrame = [];
rejectedGrid = [];

gimbalCommandTime = [];
gridAssignTime = [];
shootWorthyTime = [];
totalGridAssignTime = [];
initShootGridExec = [];

verifyGimbal = [];

for i = 1:size(data,2)
    % find where it is
    is_maxvertexN = strfind(data{i},'[Hopping_verify] MaxVertexPosN');
    is_maxvertexE = strfind(data{i},'[Hopping_verify] MaxVertexPosE');
    is_maxvertexD = strfind(data{i},'[Hopping_verify] MaxVertexPosD');
    is_minvertexN = strfind(data{i},'[Hopping_verify] MinVertexPosN');
    is_minvertexE = strfind(data{i},'[Hopping_verify] MinVertexPosE');
    is_minvertexD = strfind(data{i},'[Hopping_verify] MinVertexPosD');
    is_summary = strfind(data{i},'Grid summary');
    
    is_addedVertex = strfind(data{i},'[Hopping_verify] Added XYZ');
    
    is_gridPosImageFrame = strfind(data{i},'[Hopping_verify] gridPosImageFrame');
    is_gridPosNed = strfind(data{i},'[Hopping_verify] gridPosNed');
    
    is_rejectedGrid = strfind(data{i},'[Hopping_verify] Rejected Grid');
    is_gridSize = strfind(data{i},'[Hopping_verify] Grid size');
    
    is_gimbalCommandTime = strfind(data{i},'[Hopping_verify] Calculate gimbal commande');
    is_gridAssignTime = strfind(data{i}, '[Hopping_verify] Single grid assignment time');
    is_shootWorthyTime = strfind(data{i}, '[Hopping_verify] IsSHootWorthy');
    is_totalGridAssignTime = strfind(data{i}, '[Hopping_verify] total assign time');
    is_initShootGridExec = strfind(data{i}, '[Hopping_verify] Init Shoot Grid duration');

    is_verifyGimbal = strfind(data{i},'[Hopping_verify] [');
    
    % extract information
%     if ~isempty(is_summary)
%         is_summary
%         total_grid = data{i}(is_summary + 31:is_summar)
%         reject_grid = data{i}(is_summary + 53:end - 101)
% %         summary = vertcat(summary,str2num());
%     end
    
    if ~isempty(is_verifyGimbal)
        verifyGimbal = vertcat(verifyGimbal,str2num(data{i}(is_verifyGimbal + 17:end)));
    end
    
    if ~isempty(is_maxvertexN)
        maxvertexN = horzcat(maxvertexN,str2num(data{i}(is_maxvertexN + 31:end)));
        maxvertexN(end-9:end)=[];
    end
    if ~isempty(is_maxvertexE)
        maxvertexE = horzcat(maxvertexE,str2num(data{i}(is_maxvertexE + 31:end)));
        maxvertexE(end-9:end)=[];
    end
    if ~isempty(is_maxvertexD)
        maxvertexD = horzcat(maxvertexD,str2num(data{i}(is_maxvertexD + 31:end)));
        maxvertexD(end-9:end)=[];
    end
    if ~isempty(is_minvertexN)
        minvertexN = horzcat(minvertexN,str2num(data{i}(is_minvertexN + 31:end)));
        minvertexN(end-9:end)=[];
    end
    if ~isempty(is_minvertexE)
        minvertexE = horzcat(minvertexE,str2num(data{i}(is_minvertexE + 31:end)));
        minvertexE(end-9:end)=[];
    end
    if ~isempty(is_minvertexD)
        minvertexD = horzcat(minvertexD,str2num(data{i}(is_minvertexD + 31:end)));
        minvertexD(end-9:end)=[];
    end

    
    if ~isempty(is_addedVertex)
        addedVertexImageFrame = vertcat(addedVertexImageFrame,str2num(data{i}(is_addedVertex+29:end)));
    end 
    
    if ~isempty(is_gridPosImageFrame)
        gridPosImageFrame = vertcat(gridPosImageFrame,str2num(data{i}(is_gridPosImageFrame + 35:end)));
    end
    if ~isempty(is_gridPosNed)
        gridPosNED = vertcat(gridPosNED,str2num(data{i}(is_gridPosNed+28:end)));
    end
    
    if ~isempty(is_rejectedGrid)
        rejectedGrid = vertcat(rejectedGrid,str2num(data{i}(is_rejectedGrid+31:end)));
    end
    if ~isempty(is_gridSize)
        gridSize = str2num(data{i}(is_gridSize+29:end));
    end
    
    if ~isempty(is_gimbalCommandTime)
        gimbalCommandTime = vertcat(gimbalCommandTime,str2num(data{i}(is_gimbalCommandTime+50:end)));
    end
    if ~isempty(is_gridAssignTime)
        gridAssignTime = vertcat(gridAssignTime,str2num(data{i}(is_gridAssignTime + 47:end)));
    end
    if ~isempty(is_shootWorthyTime)
        shootWorthyTime = vertcat(shootWorthyTime, str2num(data{i}(is_shootWorthyTime + 45:end)));
    end
    if ~isempty(is_totalGridAssignTime)
        totalGridAssignTime = vertcat(totalGridAssignTime, str2num(data{i}(is_totalGridAssignTime + 37:end)));
    end
    if ~isempty(is_initShootGridExec)
        initShootGridExec = vertcat(initShootGridExec,str2num(data{i}(is_initShootGridExec+45:end)));
    end
end


%% Post processing & plotting

gridPosImageFrame = gridPosImageFrame';
gridPosNED = gridPosNED';
maxvertexNED = [maxvertexN;maxvertexE;maxvertexD];
minvertexNED = [minvertexN;minvertexE;minvertexD];
addedVertexImageFrame = addedVertexImageFrame';
addedVertexShootSurface = [addedVertexImageFrame(3,:);addedVertexImageFrame(1,:);addedVertexImageFrame(2,:)];
addedVertexNED = dcmI2ShootSurface' * addedVertexShootSurface;

maxvertexShoot = dcmI2ShootSurface * maxvertexNED;
minvertexShoot = dcmI2ShootSurface * minvertexNED;
maxvertexImage = [maxvertexShoot(2,:);maxvertexShoot(3,:);maxvertexShoot(1,:)];
minvertexImage = [minvertexShoot(2,:);minvertexShoot(3,:);minvertexShoot(1,:)];

width = gridSize(1); height = gridSize(2);
rectangle = [0 0 0 0 0 ; ...
    width/2 width/2 -width/2 -width/2 width/2 ; ...
    height/2 -height/2 -height/2 height/2 height/2];
rectangleShootSur_2 = dcmI2ShootSurface_2' * rectangle;
rectangleShootSur_1 = dcmI2ShootSurface' * rectangle;
rectangleNed_2 = [];
rectangleNed = [];
rectangleFov = rectangle;
rectangleFov(1,:) = shootDistance;

for i = 1:length(gridPosNED)
    rectangleNed{i} = rectangleShootSur_1 + gridPosNED(:,i);
end

for i = 1:length(gridPosNED)
    rectangleNed_2{i} = rectangleShootSur_2 + gridPosNED(:,i);
end

fovNed =[];
for i = 1:ftsize
    fovNed{i} = dcmI2Nav(:,:,i)' * dcmI2Gimbal(:,:,i)' * rectangleFov + posNed(:,i);
end

% vertexStartFix = 61;
% gridStartFix = 40;
% fovStartFix = 1955;

vertexStartFix = 81;
vertexEndFix = length(maxvertexD);
% vertexEndFix = 120;
gridStartFix = 65;
gridEndFix = length(gridPosNED);
% gridEndFix = 39;
fovStartFix = 65;
fovEndFix = length(fovNed);
% fovEndFix = 39;

%% plotting

% figure(1)
% clf
% plot3(maxvertexNED(2,vertexStartFix:vertexEndFix),maxvertexNED(1,vertexStartFix:vertexEndFix),-maxvertexNED(3,vertexStartFix:vertexEndFix),'r*--')
% hold on
% grid on
% plot3(minvertexNED(2,vertexStartFix:vertexEndFix),minvertexNED(1,vertexStartFix:vertexEndFix),-minvertexNED(3,vertexStartFix:vertexEndFix),'b*--')
% plot3(gridPosNED(2,gridStartFix:gridEndFix),gridPosNED(1,gridStartFix:gridEndFix),-gridPosNED(3,gridStartFix:gridEndFix),'gx-','LineWidth',0.3)
% plot3(addedVertexNED(2,:),addedVertexNED(1,:),-addedVertexNED(3,:),'kx');
% % for i = gridStartFix:gridEndFix
% %     plot3(rectangleNed_2{i}(2,:),rectangleNed_2{i}(1,:),-rectangleNed_2{i}(3,:),'k')
% % end
% for i = fovStartFix : fovEndFix
%     plot3(fovNed{i}(2,:),fovNed{i}(1,:),-fovNed{i}(3,:),'r')
% end
% axis equal
% % view(0,0)
% title('Vertex grid situation','fontsize',15)
% xlabel('E [m]','fontsize',14)
% ylabel('N [m]','fontsize',14)
% zlabel('H [m]','fontsize',14)
% legend('max vertex','min vertex','inspection order','grid position','Location','northeast')
% 
% 
% figure(2)
% clf
% plot3(gridPosImageFrame(1,:),gridPosImageFrame(3,:),-gridPosImageFrame(2,:),'ro')
% hold on
% if ~isempty(addedVertexImageFrame)
%     plot3(addedVertexImageFrame(1,:),addedVertexImageFrame(3,:),-addedVertexImageFrame(2,:),'kx')
%     hold on
% end
% plot3(maxvertexImage(1,1:end),maxvertexImage(3,1:end),-maxvertexImage(2,1:end),'bx--')
% hold on
% plot3(minvertexImage(1,1:end),minvertexImage(3,1:end),-minvertexImage(2,1:end),'rx--')
% grid on
% axis equal
% view(0,0)
% title('Vertex information','fontsize',15)
% xlabel('X [m]','FontSize',14)
% ylabel('Y [m]','fontsize',14)
% zlabel('-Z [m]','fontsize',14)
% legend('added vertex','max vertex','min vertex','Location','northeast')
% 
% % n=1;
% % for i = 1:size(data,2)
% %     if ~isnan(remem(i))
% %         tipDist(n) = str2num(data{i}(76:end-1));
% %         n = n+1;
% %     end
% % end
% 
% 
% figure(3)
% clf
% plot3(maxvertexNED(2,vertexStartFix:vertexEndFix),maxvertexNED(1,vertexStartFix:vertexEndFix),-maxvertexNED(3,vertexStartFix:vertexEndFix),'r*--')
% hold on
% grid on
% plot3(minvertexNED(2,vertexStartFix:vertexEndFix),minvertexNED(1,vertexStartFix:vertexEndFix),-minvertexNED(3,vertexStartFix:vertexEndFix),'b*--')
% plot3(gridPosNED(2,gridStartFix:gridEndFix),gridPosNED(1,gridStartFix:gridEndFix),-gridPosNED(3,gridStartFix:gridEndFix),'gx-','LineWidth',0.3)
% for i = gridStartFix:gridEndFix
%     plot3(rectangleNed{i}(2,:),rectangleNed{i}(1,:),-rectangleNed{i}(3,:),'k')
% end
% for i = fovStartFix:fovEndFix
%     plot3(fovNed{i}(2,:),fovNed{i}(1,:),-fovNed{i}(3,:),'r')
% end
% 
% plot3(ftdata.posNed_1,ftdata.posNed_0,-ftdata.posNed_2,'k:');
% plot3(ftdata.posCmdNed_1,ftdata.posCmdNed_0,-ftdata.posCmdNed_2,'r:');
% axis equal
% % view(0,0)
% title('Vertex grid situation','fontsize',15)
% xlabel('E [m]','fontsize',14)
% ylabel('N [m]','fontsize',14)
% zlabel('H [m]','fontsize',14)
% legend('max vertex','min vertex','inspection order','grid position','Location','northeast')

figure(4)
clf
plot3(maxvertexNED(2,vertexStartFix:vertexEndFix),maxvertexNED(1,vertexStartFix:vertexEndFix),-maxvertexNED(3,vertexStartFix:vertexEndFix),'r*--')
hold on
grid on
plot3(minvertexNED(2,vertexStartFix:vertexEndFix),minvertexNED(1,vertexStartFix:vertexEndFix),-minvertexNED(3,vertexStartFix:vertexEndFix),'b*--')
plot3(gridPosNED(2,gridStartFix:gridEndFix),gridPosNED(1,gridStartFix:gridEndFix),-gridPosNED(3,gridStartFix:gridEndFix),'gx-','LineWidth',0.3)
for i = gridStartFix:gridEndFix
    plot3(rectangleNed{i}(2,:),rectangleNed{i}(1,:),-rectangleNed{i}(3,:),'k')
end
for i = fovStartFix:fovEndFix
    plot3(fovNed{i}(2,:),fovNed{i}(1,:),-fovNed{i}(3,:),'r')
end

% plot3(ftdata.posNed_1,ftdata.posNed_0,-ftdata.posNed_2,'k:');
% plot3(ftdata.posCmdNed_1,ftdata.posCmdNed_0,-ftdata.posCmdNed_2,'r:');
axis equal
% view(0,0)
title('Vertex grid situation','fontsize',15)
xlabel('E [m]','fontsize',14)
ylabel('N [m]','fontsize',14)
zlabel('H [m]','fontsize',14)
legend('max vertex','min vertex','inspection order','grid position','Location','northeast')

% gimbalCommandTime
% gridAssignTime
% shootWorthyTime
% totalGridAssignTime
% initShootGridExec
