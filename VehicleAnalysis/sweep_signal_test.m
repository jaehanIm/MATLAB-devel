d2r = pi/180;
addpath 'C:\Users\dlawo\Desktop\matlab devel'

% gdLogFile = "/home/jaehan/Desktop/test flight/ref_sig_sim/gdLog_200917_135329.csv";
gdLogFile = "/home/jaehan/Desktop/test flight/ref_sig_sim/gdLog_200917_150053.csv";

[data, time] = loader(gdLogFile);
posNed = [data.posNed_0,data.posNed_1,data.posNed_2];
posNedCmd = [data.posCmdNed_0,data.posCmdNed_1,data.posCmdNed_2];
velNed = [data.velNed_0,data.velNed_1,data.velNed_2];
velNedCmd = [data.velCmdNav_0,data.velCmdNav_1,data.velCmdNav_2];
dcmI2bridge = angle2dcm(wrapToPi(data.rpy_2*d2r), zeros(size(data,1),1), zeros(size(data,1),1),'zyx');

posXyz = zeros(size(posNed));
posXyzCmd = zeros(size(posNed));
velUvw = zeros(size(posNed));
% velUvwCmd = zeros(size(posNed));
velUvwCmd = velNedCmd;
posXyzCmd = posNedCmd;
posXyz = posNed;
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

% total
angletest = 1 : length(time);
veltest = 1 : length(time);
postest = 1 : length(time);

% 배속 시뮬레이션
% angletest = 259 : 910;
% veltest = 1415 : 2067;
% postest = 2571 : 3200;

% A3 시뮬레이션
% angletest = 2553 : 3204;
% veltest = 3456 : 4108;
% postest = 4446 : 5097;


range = angletest;
figure(1)
clf
hold on
grid on
plot(time(range),data.rpy_0(range),'k')
plot(time(range),data.rpy_1(range),'b')
plot(time(range),data.rpy_2(range),'r')
plot(time(range),data.rpdCmd_0(range),'k:')
plot(time(range),data.rpdCmd_1(range),'b:')
plot(time(range),data.ySp(range),'r:')
title('angle')
size(find(islocalmax(data.rpy_0(range))),1)
legend('r','p','y')

range = veltest;
figure(2)
clf
hold on
grid on
plot(time(range),velUvw_0(range),'k')
plot(time(range),velUvw_1(range),'b')
plot(time(range),velUvw_2(range),'r')
plot(time(range),velUvwCmd_0(range),'k:')
plot(time(range),velUvwCmd_1(range),'b:')
plot(time(range),velUvwCmd_2(range),'r:')
title('vel')
size(find(islocalmax(data.velNed_1(range))),1)
legend('u','v','w')

range = postest;
figure(3)
clf
hold on
grid on
plot(time(range),posXyz_0(range),'k')
plot(time(range),posXyz_1(range),'b')
plot(time(range),posXyz_2(range),'r')
plot(time(range),posXyzCmd_0(range),'k:')
plot(time(range),posXyzCmd_1(range),'b:')
plot(time(range),posXyzCmd_2(range),'r:')
% plot(time(range),islocalmax(data.posNed_1(range))/5+mean(data.posNed_1(range)))
% plot(time(range),-islocalmin(data.posNed_1(range))/5+mean(data.posNed_1(range)))
% plot(time(range),data.rpy_0(range)/10+mean(data.posNed_1(range)))
% plot(time(range),data.velNed_1(range)+mean(data.posNed_1(range)))
title('pos')
size(find(islocalmax(data.posNed_1(range))),1)
legend('x','y','z')

% figure(4)
% clf
% hold on
% grid on
% plot(time(range),velUvwCmd_1(range)-velUvw_1(range))
% plot(time(range),velUvwCmd_1(range))


figure(5)
clf
hold on
grid on
title('response time plot')
plot(time,data.fcMcMode)

%% Reference sweep signal

% figure(5)
% clf
% c1 = 4; c2 = 0.01866;
% wmin = 0.4*2*pi;
% wmax = 10*2*pi;
% t = 0:0.02:13;
% T = 13;
% res = c2*(wmax-wmin)*T/c1;
% omega = wmin + (exp(t/T*c1)-1)*c2*(wmax-wmin);
% theta = wmin*t + c2*(wmax-wmin)*(T/c1*exp(c1/T*t)-t) - res;
% plot(t,sin(theta))
% plot(t,omega/2/pi)