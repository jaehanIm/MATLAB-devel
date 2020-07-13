%% Parameter setting
TOPos = [10 20 0];        % T/O GPS position (including T/O field elevation)0
currPos_rel = [0 0 2]; % Relative GPS NED to T/O position (after auto T/O)
currPos = [];           % GPS position of drone
wtPos = [100 100];      % Absolute GPS position of WT
wtHeight = 70;          % WT nose height from 'ground'
wtElev = 10;            % WT installation field elevation
wtAngle = 60;           % yaw angle in degrees
bladeLen = 30;          % WT blade length
noFlyWidth = 5;

% initialization
isNear_ = false;         % Escape maneuver flag
isRound_ = false;        % Round trajectory flag
horDist = 0;           % Horizontal distance between drone and WT
vert_dist = 0;          % Vertical distance between drone and WT
currPos = TOPos + currPos_rel;
safeDist = bladeLen;    % Safe distance - 일단 blade length 만큼 회피
devDistHor = 0;            % Deviation distance from center line
devDistLon = 0;
wtDirec = -atan2(wtPos(1)-currPos(1),wtPos(2)-currPos(2)) + pi/2;
wtDirec = mod(wtDirec,2*pi);
coreAngle = 0;

d2r = pi/180;
wtAngle = wtAngle * d2r;

wp = [];
pathType = 0;
isFinal_ = false;
isNoFly_ = false;

%% Type discrimination
horDist = norm(currPos(1:2)-wtPos);
thres_dist = sqrt((horDist)^2-(safeDist)^2);
interPos = wtPos + [safeDist*cos(wtAngle) safeDist*sin(wtAngle)];
inter_dist = norm(interPos - currPos(1:2));
interDirec = -atan2(interPos(1)-currPos(1),interPos(2)-currPos(2)) + pi/2;
interDirec = mod(interDirec,2*pi);
coreAngle = acos(safeDist/horDist);

% if escape maneuver is required
if horDist < safeDist
    isNear_ = true;
    % Escape wp type
    [isToFrom_, devDistLon] = isToFrom(currPos(1:2),wtPos,wtAngle);    
    [isLR_, devDistHor] = isLR(currPos(1:2),wtPos,wtAngle);
    isNoFly_ = isNoFly(currPos(1:2),wtPos,wtAngle,noFlyWidth);
else
    isNear_ = false;
end 

% if round trajectory required
if thres_dist >= inter_dist
    isRound_ = false;
else
    isRound_ = true;
end

%% WP generation (phase 1)
if isNear_ % Inside safety zone
    if isToFrom_ == "front"
        pathType = 1;
        isFinal_ = 1;
        wp = vertcat(wp, [interPos, wtHeight, pathType, isFinal_]);
    else % Located at the back of WT -> Most complex situation
        if isLR_ == "left"
            direction = wtAngle + pi/2;
        elseif isLR_ == "right"
            direction = wtAngle - pi/2;
        end
        residue = sqrt(bladeLen^2 - devDistLon^2);
        midPos = currPos(1:2) + (residue - devDistHor) .* [cos(direction) sin(direction)];
        pathType = 1;
        isFinal_ = 0;
        wp = vertcat(wp, [midPos, wtHeight, pathType, isFinal_]);
    end
else % Out of safety zone
    if isRound_ % Requires Circular Trajectory
        [~,type] = anglerel(wtDirec,interDirec);
        if type == "right"
            direction = wtDirec + (coreAngle - pi/2);
            midPos = thres_dist .* [cos(direction) sin(direction)] + currPos(1:2);
        elseif type == "left"
            direction = wtDirec - (coreAngle - pi/2);
            midPos = thres_dist .* [cos(direction) sin(direction)] + currPos(1:2);
        end
        pathType = 1;
        isFinal_ = false;
        wp = vertcat(wp, [midPos, wtHeight, pathType, isFinal_]);
    else % Simple Straight Trajectory
        pathType = 1;
        isFinal_ = true;
        wp = vertcat(wp, [interPos, wtHeight, pathType,isFinal_]);
    end
end
%% WP generation (phase 2)
if isFinal_ == 1           % yaw angle in degrees
    disp("WP generation complete")
else 
    pathType = 2;
    isFinal_ = true;
    wp = vertcat(wp, [interPos, wtHeight, pathType, isFinal_]);
    disp("WP generation complete")
end

%% Trajectory generation
if isNoFly_ % Inappropriate drone location -> Remove trajectory
    wp = [];
    disp('Drone is located in NoFlyZone')
    disp('=== Terminating generated information ===')
else
    trajectory = [];
    trajectory = vertcat(trajectory,currPos);
    for i = 1 : size(wp,1)
        if wp(i,4) == 1
            trajectory = vertcat(trajectory,wp(i,1:3));
        elseif wp(i,4) == 2
            mid_direc = -atan2(midPos(1)-wtPos(1),midPos(2)-wtPos(2)) + pi/2;
            inter_direc = wtAngle;
            tempFlag_ = isLR(interPos,wtPos,mid_direc);
            
            if tempFlag_ == "right"
                [temp_x, temp_y] = draw_circle(wtPos,safeDist,'r',[mid_direc inter_direc]);
            elseif tempFlag_ == "left"
                [temp_x, temp_y] = draw_circle(wtPos,safeDist,'r',[mid_direc inter_direc]);
            end
            
            trajectory = vertcat(trajectory,[temp_x', temp_y', ones(length(temp_x),1).*wtHeight]);
        end
    end
end

if isNoFly_ == false
%% Plot
figure(1)
% clf
hold on
plot(currPos(1),currPos(2),'*')
plot(wtPos(1),wtPos(2),'Xr','MarkerSize',10,'LineWidth',2)
% plot([wtPos(1) currPos(1)],[wtPos(2) currPos(2)],'b:')
draw_circle(wtPos,safeDist,'r');
% plot([currPos(1) interPos(1)],[currPos(2) interPos(2)],'g:')
plot(interPos(1),interPos(2),'xm','MarkerSize',10','LineWidth',3)
NoFlyBox = [noFlyWidth -noFlyWidth -noFlyWidth noFlyWidth, noFlyWidth; ...
    bladeLen, bladeLen, -bladeLen, -bladeLen, bladeLen];
NoFlyBox = [cos(wtAngle), -sin(wtAngle) ; sin(wtAngle), cos(wtAngle)]*NoFlyBox;
NoFlyBox = NoFlyBox + [wtPos(1) ; wtPos(2)];
plot(NoFlyBox(1,:),NoFlyBox(2,:),'k--')
grid on
plot_arrow(wtPos(1), wtPos(2), wtPos(1)-10*cos(wtAngle+pi), wtPos(2)-10*sin(wtAngle+pi),'linewidth',2,'color',[0.5 0.5 0.5]);

for i = 1:size(wp)
    plot(wp(i,1),wp(i,2),'bo','LineWidth',3)
end
plot(wp(1,1),wp(1,2),'ro','LineWidth',3)

plot(trajectory(:,1), trajectory(:,2), 'k')

title('WT Path Finding','fontsize',14)
xlabel('East [m]','fontsize',14)
ylabel('North [m]','FontSize',14)

%% 3d PLOT
% figure(2)
% % clf
% hold on
% plot3(currPos(1),currPos(2),currPos(3),'*')
% plot3(wtPos(1),wtPos(2),wtHeight,'Xr','MarkerSize',10,'LineWidth',2)
% % plot3([wtPos(1) currPos(1)],[wtPos(2) currPos(2)],[wtHeight currPos(3)],'b:')
% draw_circle3([wtPos wtHeight],safeDist,'r');
% % plot3([currPos(1) interPos(1)],[currPos(2) interPos(2)],[currPos(3) wtHeight],'g:')
% plot3(interPos(1),interPos(2),wtHeight,'xm','MarkerSize',10','LineWidth',3)
% NoFlyBox = [noFlyWidth -noFlyWidth -noFlyWidth noFlyWidth, noFlyWidth; ...
%     bladeLen, bladeLen, -bladeLen, -bladeLen, bladeLen];
% NoFlyBox = [cos(wtAngle), -sin(wtAngle) ; sin(wtAngle), cos(wtAngle)]*NoFlyBox;
% NoFlyBox = NoFlyBox + [wtPos(1) ; wtPos(2)];
% plot3(NoFlyBox(1,:),NoFlyBox(2,:),ones(1,5)*wtHeight,'k--')
% grid on
% 
% for i = 1:size(wp)
%     plot3(wp(i,1),wp(i,2),wp(i,3),'bo','LineWidth',3)
% end
% plot3(wp(1,1),wp(1,2),wp(1,3),'ro','LineWidth',3)
% 
% plot3(trajectory(:,1), trajectory(:,2), trajectory(:,3), 'k')
% 
% title('WT Path Finding','fontsize',14)
% xlabel('East [m]','fontsize',14)
% ylabel('North [m]','FontSize',14)
end
