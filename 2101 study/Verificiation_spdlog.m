clear all
% close all

path = '/home/jaehan/log/sample6_pitch90_blade0 whole_vertexaddition update/spdlog_210106_165604.txt';
% gdpath = '/home/jaehan/log/210106_120246/gdLogCsv/gdLog_210106_120246.csv';

fileID = fopen(path);
% [ftdata,time] = loader(gdpath);

n = 1;
while ~feof(fileID)
    data{n} = fgetl(fileID);
    n = n+1;
end

d2r = pi/180;
dcmI2ShootSurface = angle2dcm(0,0,0,'zyx');

maxvertexN = [];
maxvertexE = [];
maxvertexD = [];
minvertexN = [];
minvertexE = [];
minvertexD = [];
gridPosImageFrame = [];
gridPosNED = [];
addedVertexImageFrame = [];
rejectedGrid = [];

for i = 1:size(data,2)
    % find where it is
    is_maxvertexN = strfind(data{i},'[Hopping_verify] MaxVertexPosN');
    is_maxvertexE = strfind(data{i},'[Hopping_verify] MaxVertexPosE');
    is_maxvertexD = strfind(data{i},'[Hopping_verify] MaxVertexPosD');
    is_minvertexN = strfind(data{i},'[Hopping_verify] MinVertexPosN');
    is_minvertexE = strfind(data{i},'[Hopping_verify] MinVertexPosE');
    is_minvertexD = strfind(data{i},'[Hopping_verify] MinVertexPosD');
    
    is_addedVertex = strfind(data{i},'[Hopping_verify] Added XYZ');
    
    is_gridPosImageFrame = strfind(data{i},'[Hopping_verify] gridPosImageFrame');
    is_gridPosNed = strfind(data{i},'[Hopping_verify] gridPosNed');
    
    is_rejectedGrid = strfind(data{i},'[Hopping_verify] Rejected Grid');
    is_gridSize = strfind(data{i},'[Hopping_verify] Grid size');
    
    % extract information
    if ~isempty(is_maxvertexN)
        maxvertexN = horzcat(maxvertexN,str2num(data{i}(is_maxvertexN + 31:end)));
    end
    if ~isempty(is_maxvertexE)
        maxvertexE = horzcat(maxvertexE,str2num(data{i}(is_maxvertexE + 31:end)));
    end
    if ~isempty(is_maxvertexD)
        maxvertexD = horzcat(maxvertexD,str2num(data{i}(is_maxvertexD + 31:end)));
    end
    if ~isempty(is_minvertexN)
        minvertexN = horzcat(minvertexN,str2num(data{i}(is_minvertexN + 31:end)));
    end
    if ~isempty(is_minvertexE)
        minvertexE = horzcat(minvertexE,str2num(data{i}(is_minvertexE + 31:end)));
    end
    if ~isempty(is_minvertexD)
        minvertexD = horzcat(minvertexD,str2num(data{i}(is_minvertexD + 31:end)));
    end
    
    if ~isempty(is_addedVertex) %% TODO
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
end

gridPosImageFrame = gridPosImageFrame';
gridPosNED = gridPosNED';
maxvertexNED = [maxvertexN;maxvertexE;maxvertexD];
minvertexNED = [minvertexN;minvertexE;minvertexD];
addedVertexImageFrame = addedVertexImageFrame';

maxvertexShoot = dcmI2ShootSurface * maxvertexNED;
minvertexShoot = dcmI2ShootSurface * minvertexNED;
maxvertexImage = [maxvertexShoot(2,:);maxvertexShoot(3,:);maxvertexShoot(1,:)];
minvertexImage = [minvertexShoot(2,:);minvertexShoot(3,:);minvertexShoot(1,:)];

width = gridSize(1); height = gridSize(2);
rectangle = [0 0 0 0 0 ; ...
    width/2 width/2 -width/2 -width/2 width/2 ; ...
    height/2 -height/2 -height/2 height/2 height/2];
rectangle = dcmI2ShootSurface' * rectangle;
rectangleNed = [];
for i = 1:length(gridPosNED)
    rectangleNed{i} = rectangle + gridPosNED(:,i);
end

figure(1)
clf
plot3(maxvertexNED(2,:),maxvertexNED(1,:),-maxvertexNED(3,:),'r*--')
hold on
grid on
plot3(minvertexNED(2,:),minvertexNED(1,:),-minvertexNED(3,:),'b*--')
plot3(gridPosNED(2,:),gridPosNED(1,:),-gridPosNED(3,:),'gx-','LineWidth',0.3)
for i = 1:length(gridPosNED)
    plot3(rectangleNed{i}(2,:),rectangleNed{i}(1,:),-rectangleNed{i}(3,:),'k')
end
axis equal
view(0,0)
title('Vertex grid situation','fontsize',15)
xlabel('E [m]','fontsize',14)
ylabel('N [m]','fontsize',14)
zlabel('H [m]','fontsize',14)
legend('max vertex','min vertex','inspection order','grid position')

rejectedGrid+1;
figure(2)
clf
% plot3(gridPosImageFrame(1,:),gridPosImageFrame(3,:),-gridPosImageFrame(2,:),'ro')
% hold on
if ~isempty(addedVertexImageFrame)
    plot3(addedVertexImageFrame(1,:),addedVertexImageFrame(3,:),-addedVertexImageFrame(2,:),'kx')
    hold on
end
plot3(maxvertexImage(1,:),maxvertexImage(3,:),-maxvertexImage(2,:),'bx--')
hold on
plot3(minvertexImage(1,:),minvertexImage(3,:),-minvertexImage(2,:),'rx--')
grid on
axis equal
view(0,0)
title('Vertex information','fontsize',15)
xlabel('X [m]','FontSize',14)
ylabel('Y [m]','fontsize',14)
zlabel('-Z [m]','fontsize',14)
legend('added vertex','max vertex','min vertex')

% n=1;
% for i = 1:size(data,2)
%     if ~isnan(remem(i))
%         tipDist(n) = str2num(data{i}(76:end-1));
%         n = n+1;
%     end
% end
% 
