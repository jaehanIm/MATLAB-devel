clear all
% close all

directory = "/home/jaehan/Desktop/white_creek_additional/";
% directory = "/home/jaehan/Desktop/white_creek_sample/";
list = dir(directory);
d2r = pi/180; r2d = 180/pi;
R_data = [];

for filenum = 3 : size(list,1)
% for filenum = 89

% filenum = 21;
data = readtable(directory+list(filenum).name);
N = data.n; E = data.e; D = data.alt; yaw = data.bodyYaw * d2r; gimbp = data.gimbalPitch * d2r;
blnum = data.bladePosition; blside = data.bladeSide; direction = data.direction;
curPosNed = [N,E,D];

startval = 0; temp = blnum(1);
while temp == 255 && startval ~= size(data,1)
    startval = startval + 1;
    temp = blnum(startval);
end
if startval == size(data ,1)
    
end
if startval == 0
    startval = 1;
end
wtNosePosNed = [N(startval),E(startval),D(startval)];
wtTailInitYaw = yaw(startval);

if blnum(startval) == 0
    mission = "6-1";
elseif blnum(startval) == 1
    mission = "6-2";
else
    mission = "Unidentified";
end

%% job categorization
job = cell(size(data,1),1);
jobnum = zeros(size(data,1),1);
for i = 1:size(data,1)    
    if length(direction{i}) == 9 && blnum(i) == 2 && (blside{i} == "TE" || length(blside{i}) == 9)
        job{i} = 'wp'; jobnum(i) = 2;
    elseif length(direction{i}) == 9 && blnum(i) ~= 255
        job{i} = 'ct'; jobnum(i) = 1;
    elseif blnum(i) ~= 255
        job{i} = 'nm'; jobnum(i) = 0;
    elseif blnum(i) == 255
        if i > 1
            job{i} = job{i-1}; jobnum(i) = -1;
            blnum(i) = blnum(i-1);
        else
            job{i} = 'tt'; jobnum(i) = -1;
        end
    end
end

%% R calculation
R = zeros(size(data,1),1);
dcmI2WtTail = angle2dcm(wrapToPi(wtTailInitYaw), 0, 0,'zyx');

for num = 1 : size(data,1)
    % basic R calculation
    curPosNedWtNose = curPosNed(num,:) - wtNosePosNed;
    if blnum(num) == 1
        R(num) = abs(curPosNedWtNose(3));
    else        
        curPosXyzWtNose = dcmI2WtTail * curPosNedWtNose';
        R(num) = curPosXyzWtNose(3)*sin(30*d2r) + abs(curPosXyzWtNose(2))*cos(30*d2r);
    end
    
    % R value adjustment
    inspectc = 7;
    if job{num} == 'wp'
        yawDevWTAngle = abs(wrapToPi(wtTailInitYaw+pi - yaw(num)));
        R(num) = (abs(curPosXyzWtNose(2))-abs(curPosXyzWtNose(1)) * tan(abs(yawDevWTAngle)))/cos(30*d2r);
    elseif job{num} == 'ct'
        if blnum(num) == 1
            correctedDist = inspectc * tan(gimbp(num));
            R(num) = R(num) + correctedDist;
        else
            gimbalPitchDevAngle = 60*d2r + gimbp(num);
            correctedDist = inspectc * tan(gimbalPitchDevAngle);
            R(num) = R(num) - correctedDist;
        end
    end
    
    % Negative value protection
    if R(num) < 0
        R(num) = 0;
    end
    
    % Test value removal
    if job{num} == 'tt'
        R(num) = 0;
    end
    
end

note = [];
if mission == "Unidentified"
    note = ' unidentified mission';
end

R_data{filenum-2} = R;

%% Write

data.R = R;
writetable(data,directory+list(filenum).name);

%% plotting
% figure()
% hold on
% plot(R,'.')
% plot(jobnum*10)
% % ylim([0 50]);
% plot(R,'b:')
% legend('R','job type')
% xlim([0 size(data,1)])
% if mission == "6-2"
%     text(90,45,"Mission 6-2")
%     text(35, 15, "bl #1")
%     text(170, 15, "bl #1")
%     text(260, 15, "bl #2")
%     text(90,40,list(filenum).name);
%     text(90,35,num2str(filenum));
% elseif mission == "6-1"
%     text(90,45,"Mission 6-1")
%     text(35, 15, "bl #0")
%     text(170, 15, "bl #1")
%     text(260, 15, "bl #1")
%     text(90,40,list(filenum).name);
%     text(90,35,num2str(filenum));
% else
%     text(90,45,"Unknown")
%     text(90,40,list(filenum).name);
%     text(90,35,num2str(filenum));
% end


%% Check writing 
data2 = readtable(directory+list(filenum).name);
note2 = []; note3 = [];
if any((abs(R - data2.R)<1e-5)-1)
    note2 = ' / bad writing detected'; 
end
if length(R) < 280
    note3 = ' / mission abort expected';
end

disp(['Processed file: ',list(filenum).name,' (',num2str(filenum),'/',num2str(size(list,1)),')',note,note2,note3]);

end