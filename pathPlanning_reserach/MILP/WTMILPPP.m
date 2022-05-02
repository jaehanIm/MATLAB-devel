d2r = pi/180;
r2d = 1/d2r;

%% WT param.

% NED coord.
wtPos = [0 0 0]';
wtDir = 0; %deg (Looking into the blade)
wtHeight = 562;
bladeL = 255;
noseL = 50;
homePos = [-30 10 0]';
wtRot = 0;
inspectionDist = 20;
rootDist = 30;
vnum = 1;

stlAddr = '/home/jaehan/Desktop/WT.stl';

%% Inspection point generation 

wtDir = wtDir * d2r;
wtRot = wtRot * d2r;
dcmI2BlPlane = angle2dcm(wtDir,0,wtRot,'zyx');

dcmBlPlane2Bl = [];
dcmBlPlane2Bl{1} = angle2dcm(0,0,120*d2r,'zyx');
dcmBlPlane2Bl{2} = angle2dcm(0,0,0*d2r,'zyx');
dcmBlPlane2Bl{3} = angle2dcm(0,0,-120*d2r,'zyx');

wtCore = wtPos + [0 0 -wtHeight]';
nosePos = wtCore + dcmI2BlPlane' * [-noseL, 0, 0]';

rootWp = zeros(3,4,3); % bl#, fix#, posNed
tipWp = zeros(3,4,3);
centerWp = zeros(3,4,3);
for blNum = 1:3
    tempRoot = [];
    tempTip = [];
    tempCenter = [];
    for fixNum = 0:3
        tempRoot(fixNum+1,:) = [inspectionDist * cos(pi/2*fixNum), inspectionDist * sin(pi/2*fixNum), -rootDist];
        tempTip(fixNum+1,:) = [inspectionDist * cos(pi/2*fixNum), inspectionDist * sin(pi/2*fixNum), -bladeL];
        tempCenter(fixNum+1,:) = [inspectionDist * cos(pi/2*fixNum), inspectionDist * sin(pi/2*fixNum), -(bladeL+rootDist)/2];
    end
    rootWp(blNum,:,:) = (wtCore + dcmI2BlPlane' * dcmBlPlane2Bl{blNum}' * tempRoot')';
    tipWp(blNum,:,:) = (wtCore + dcmI2BlPlane' * dcmBlPlane2Bl{blNum}' * tempTip')';
    centerWp(blNum,:,:) = (wtCore + dcmI2BlPlane' * dcmBlPlane2Bl{blNum}' * tempCenter')';
end

totalWp(1,:,:,:) = rootWp;
totalWp(2,:,:,:) = tipWp;
totalWp(3,:,:,:) = centerWp;

figure(1)
clf
hold on
grid on
for j = 1:3
tempRoot = reshape(rootWp(j,:,:),4,3);
tempTip = reshape(tipWp(j,:,:),4,3);
tempCenter = reshape(centerWp(j,:,:),4,3);
for i = 1:4
plot3(tempRoot(i,2),tempRoot(i,1),-tempRoot(i,3),'k*')
plot3(tempTip(i,2),tempTip(i,1),-tempTip(i,3),'k*')
% plot3(tempCenter(i,2),tempCenter(i,1),-tempCenter(i,3),'k*')
end
end

plot3(wtCore(2),wtCore(1),-wtCore(3),'rx')
plot3(nosePos(2),nosePos(1),-nosePos(3),'o','LineWidth',3)
plot3(homePos(2),homePos(1),homePos(3),'o','MarkerSize',5,'LineWidth',3)
view(35,40)
axis equal

%% Network formulation

% A -> connectivity matrix
% C -> cost matrix
N = 24 + 1; % add nose
A = zeros(N,N);
C = zeros(N,N);
IMR = [];

% root to tip connector
for blNum = 1:3
    base = (blNum-1)*8 + 1;
    self = [0 1 0 1;1 0 1 0];
    self = vertcat(self,self);
    rootTip = eye(4);
    region = [self,rootTip;rootTip,self];
    A(base+1:base+size(region,1),base+1:base+size(region,1)) = region;
end

% inter blade connector
region = zeros(8,8);
region12 = region;
region12(1:4,1:4) = ...
    [0 0 0 0;...
    0 0 0 0;...
    0 0 1 0;...
    0 1 0 0];
region13 = region;
region13(1:4,1:4) = ...
    [0 0 0 0;...
    0 0 0 0;...
    0 0 1 0;...
    0 0 0 0];
region23 = region12;
A(2:9,10:25) = [region12,region13];
A(10:17,18:25) = region23;
A(10:25,2:9) = [region12,region13]';
A(18:25,10:17) = region23';

% nose blade connector
temp = zeros(1,N);
for blNum = 1:3
    base = (blNum-1)*8 + 1;
    temp(base+3) = 1;
end
A(1,:) = temp;
A(:,1) = temp';

%% Cost matrix formulation
% trivial ones
for i = 1:N
    for j = 1:N
        if A(i,j) == 1 && i~=j
            initNod = nodeInterpreter(i);
            targNod = nodeInterpreter(j);
            if initNod.home
                startNodePos = nosePos;
            else
                startNodePos = reshape(totalWp(initNod.rt,initNod.blNum,initNod.fixNum,:),3,1);
            end
            if targNod.home
                endNodePos = nosePos;
            else
                endNodePos = reshape(totalWp(targNod.rt,targNod.blNum,targNod.fixNum,:),3,1);
            end
            C(i,j) = norm(endNodePos - startNodePos);
        end
    end
end

wtGraph = graph(C);

figure(3)
clf
plot(wtGraph)

% complet-ification
implicitRoute = [];
for i = 1:N
    for j = i:N
        if A(i,j) ~= 1 && i~=j
            [implRoute,d]=shortestpath(wtGraph,i,j,'Method','positive');
            implicitRoute{i,j} = implRoute;
            implicitRoute{j,i} = fliplr(implRoute);
            C(i,j) = d;
            C(j,i) = d;
        end
    end
end

%% Solve!
map.N = N;
map.capacity = bladeL * 12 / vnum;
% map.capacity = 1e10;
map.C = C;
map.vnum = vnum;
map.oblig = ...
    [4,8;3,7;2,6;5,9;...
    12,16;11,15;10,14;13,17;...
    20,24;21,25;18,22;19,23];
solve_complete = false;
tic
while ~solve_complete
    [routeResult,score] = TSP_solver(map);
    if score == -1
        map.capacity = map.capacity * 1.1;
    else
        solve_complete = true;
    end
end
toc

% Home edition
homeAddNode = [];
for i = 1:size(routeResult,1)
    if routeResult(i,end) ~= 0
        homeAddNode(i) = 1;
    else
        zeroRoute = find(routeResult(i,:)==0);
        zeroRouteIdx = zeroRoute(1);
        routeResult(i,zeroRouteIdx) = 1;
        homeAddNode(i) = 0;
    end
end
routeResult = horzcat(routeResult,homeAddNode');

% Route extraction
Nr = size(routeResult,2);
routePos = zeros(map.vnum,3,Nr);
for v = 1:map.vnum
    for i = 1:size(routeResult,2)
        routeNod = nodeInterpreter(routeResult(v,i));
        if routeNod.home
            routePos(v,:,i) = nosePos;
        else
            routePos(v,:,i) = reshape(totalWp(routeNod.rt,routeNod.blNum,routeNod.fixNum,:),3,1);
        end
    end
end

% Route decipher
routeResultDec = [];
for v = 1:map.vnum
    routeResultDec{v} = routeResult(v,1);
    addRoute = [];
    for i = 1:Nr-1
        ni = routeResult(v,i);
        nf = routeResult(v,i+1);
        if nf == 0
            addRoute = 0;
        else
            if A(ni,nf) == 0
                hiddenRoute = implicitRoute{ni,nf};
                hiddenRoute(1) = [];
                addRoute = hiddenRoute;
            else
                addRoute = nf;
            end
        end
        routeResultDec{v}(end+1:end+length(addRoute)) = addRoute;
    end
end

% Deciphered route extraction
routePosDec = [];
for v = 1:map.vnum
    Nd = size(routeResultDec{v},2);
    for i = 1:Nd
        routeNod = nodeInterpreter(routeResultDec{v}(i));
        if routeNod.home
            routePosDec{v}(:,i) = nosePos;
        else
            routePosDec{v}(:,i) = reshape(totalWp(routeNod.rt,routeNod.blNum,routeNod.fixNum,:),3,1);
        end
    end
end

%% Plot
figure(11)
clf
hold on
grid on
for j = 1:3
tempRoot = reshape(rootWp(j,:,:),4,3);
tempTip = reshape(tipWp(j,:,:),4,3);
for i = 1:4
plot3(tempRoot(i,2),tempRoot(i,1),-tempRoot(i,3),'k*')
plot3(tempTip(i,2),tempTip(i,1),-tempTip(i,3),'k*')
end
end
plot3(wtCore(2),wtCore(1),-wtCore(3),'rx')
plot3(nosePos(2),nosePos(1),-nosePos(3),'mo','MarkerSize',3,'LineWidth',7)
plot3(homePos(2),homePos(1),homePos(3),'o','MarkerSize',1,'LineWidth',3)
view(35,40)

axis equal

% route plot
for v = 1:map.vnum
    plot3(reshape(routePos(v,2,:),1,Nr),reshape(routePos(v,1,:),1,Nr),-reshape(routePos(v,3,:),1,Nr),'*--','LineWidth',1)
end
title('Raw route result')

figure(12)
clf
hold on
grid on
for j = 1:3
tempRoot = reshape(rootWp(j,:,:),4,3);
tempTip = reshape(tipWp(j,:,:),4,3);
for i = 1:4
plot3(tempRoot(i,2),tempRoot(i,1),-tempRoot(i,3),'k*')
plot3(tempTip(i,2),tempTip(i,1),-tempTip(i,3),'k*')
end
end
plot3(wtCore(2),wtCore(1),-wtCore(3),'rx')
plot3(nosePos(2),nosePos(1),-nosePos(3),'mo','MarkerSize',3,'LineWidth',7)
plot3(homePos(2),homePos(1),homePos(3),'o','MarkerSize',1,'LineWidth',3)
view(35,40)
axis equal

% deciphered route plot
for v = 1:map.vnum
    plot3(routePosDec{v}(2,:),routePosDec{v}(1,:),-routePosDec{v}(3,:),'*--','LineWidth',1)
end
title('Final result')