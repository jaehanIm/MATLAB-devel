gdLogFile = "/home/jaehan/Desktop/test flight/Vehicle_Analysis/KH/200904_141145/gdLog_200904_141145.csv";
% gdLogFile = "/home/jaehan/Desktop/test flight/Vehicle_Analysis/SB/200904_144201/gdLog_200904_144201.csv";

[data, data_time] = loader(gdLogFile);
data_time = seconds(data_time);
d2r = pi/180;

% Sweep signal param.
c1 = 4; c2 = 0.01866;
wmin = 0.4*2*pi;
wmax = 10*2*pi;
T = 20;

%% Coordination adjustment
posNed = [data.posNed_0,data.posNed_1,data.posNed_2];
posNedCmd = [data.posCmdNed_0,data.posCmdNed_1,data.posCmdNed_2];
velNed = [data.velNed_0,data.velNed_1,data.velNed_2];
velNedCmd = [data.velCmdNav_0,data.velCmdNav_1,data.velCmdNav_2];
dcmI2bridge = angle2dcm(wrapToPi(data.rpy_2*d2r), zeros(size(data,1),1), zeros(size(data,1),1),'zyx');

posXyz = zeros(size(posNed));
posXyzCmd = zeros(size(posNed));
velUvw = zeros(size(posNed));
velUvwCmd = velNedCmd;
% posXyzCmd = posNedCmd;
% posXyz = posNed;
for i = 1:size(posNed,1)
    posXyz(i,:) = dcmI2bridge(:,:,i) * posNed(i,:)';
    posXyzCmd(i,:) = dcmI2bridge(:,:,i) * posNedCmd(i,:)';
    velUvw(i,:) = dcmI2bridge(:,:,i) * velNed(i,:)';
%     velUvwCmd(i,:) = dcmI2bridge(:,:,i) * velNedCmd(i,:)';
end

posXyz_0 = posXyz(:,1);
posXyz_1 = posXyz(:,2);
posXyz_2 = posXyz(:,3);
posXyzCmd_0 = posXyzCmd(:,1);
posXyzCmd_1 = posXyzCmd(:,2);
posXyzCmd_2 = posXyzCmd(:,3);
velUvw_0 = velUvw(:,1);
velUvw_1 = velUvw(:,2);
velUvw_2 = velUvw(:,3);
velUvwCmd_0 = velUvwCmd(:,1);
velUvwCmd_1 = velUvwCmd(:,2);
velUvwCmd_2 = velUvwCmd(:,3);

%% Data parsing 
testStartFlag = [];
testFinishFlag = [];
missionStartFlag = [];
missionFinishFlag = [];
fcMode = data.fcMcMode;
missionType = data.missionType;

ns = 1; nf = 1; nsm = 1; nfm = 1;
for i = 2:length(fcMode)
    if fcMode(i) ~= fcMode(i-1) && fcMode(i) == 2
        testStartFlag(ns) = i; ns = ns+1;
    elseif fcMode(i) ~= fcMode(i-1) && fcMode(i-1) == 2
        testFinishFlag(nf) = i; nf = nf+1;
    end
    if fcMode(i-1) == 0 && fcMode(i) == 1
        missionStartFlag(nsm) = i; nsm = nsm + 1;
    elseif fcMode(i-1) == 1 && fcMode(i) == 0
        missionFinishFlag(nfm) = i; nfm = nfm + 1;
    end
end

%% Data Selection
% Mission Type
% 4-1: 2 / 4-2: 3 / 4-3: 4 / 5-1: 0 / 5-2: 1 / 6-1: 8 / 6-2: ?

n = 16;

% Cmd = data.rpdCmd_0(testStartFlag(n):testFinishFlag(n));
Cmd = data.ySp(testStartFlag(n):testFinishFlag(n));
SigStartFlag = find(Cmd);
SigFinishFlag = SigStartFlag(end)+testStartFlag(n)-2;
SigStartFlag = SigStartFlag(1)+testStartFlag(n);

Response = data.rpy_0(SigStartFlag:SigFinishFlag);
Cmd = data.rpdCmd_0(SigStartFlag:SigFinishFlag);
% Response = posXyz_1(SigStartFlag:SigFinishFlag);
% Cmd = posXyzCmd_1(SigStartFlag:SigFinishFlag);
time = data_time(SigStartFlag:SigFinishFlag);

res = c2*(wmax-wmin)*T/c1;
% omega = wmin + (exp((time-time(1))/T*c1)-1)*c2*(wmax-wmin);
freq = omega/2/pi;
% theta = wmin*t + c2*(wmax-wmin)*(T/c1*exp(c1/T*t)-t) - res;


%% Plotting

range = SigStartFlag:SigFinishFlag;

figure(1)
clf
hold on
grid on
plot(time,data.rpy_0(range),'k')
plot(time,data.rpy_1(range),'b')
plot(time,data.rpy_2(range),'r')
plot(time,data.rpdCmd_0(range),'k:')
plot(time,data.rpdCmd_1(range),'b:')
plot(time,data.ySp(range),'r:')
title('angle')
size(find(islocalmax(data.rpy_0(range))),1)
legend('r','p','y')

figure(2)
clf
hold on
grid on
plot(time,velUvw_0(range),'k')
plot(time,velUvw_1(range),'b')
plot(time,velUvw_2(range),'r')
plot(time,velUvwCmd_0(range),'k:')
plot(time,velUvwCmd_1(range),'b:')
plot(time,velUvwCmd_2(range),'r:')
title('vel')
size(find(islocalmax(data.velNed_1(range))),1)
legend('u','v','w')

figure(3)
clf
hold on
grid on
plot(time,posXyz_0(range),'k')
plot(time,posXyz_1(range),'b')
plot(time,posXyz_2(range),'r')
plot(time,posXyzCmd_0(range),'k:')
plot(time,posXyzCmd_1(range),'b:')
plot(time,posXyzCmd_2(range),'r:')
% plot(time,islocalmax(data.posNed_1(range))/5+mean(data.posNed_1(range)))
% plot(time,-islocalmin(data.posNed_1(range))/5+mean(data.posNed_1(range)))
% plot(time,data.rpy_0(range)/10+mean(data.posNed_1(range)))
% plot(time,data.velNed_1(range)+mean(data.posNed_1(range)))
title('pos')
size(find(islocalmax(data.posNed_1(range))),1)
legend('x','y','z')

% figure(4)
% clf
% hold on
% grid on
% title('response time plot')
% plot(time,data.fcMcMode(range))