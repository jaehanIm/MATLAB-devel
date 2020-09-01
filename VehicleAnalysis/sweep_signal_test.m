addpath 'C:\Users\dlawo\Desktop\matlab devel'

% gdLogFile = "/home/jaehan/Desktop/test flight/ref_sig_sim/gdLog_200813_121452.csv"; % 배속 시뮬레이션 데이터
gdLogFile = "/home/jaehan/Desktop/test flight/ref_sig_sim/gdLog_200812_152705.csv"; % A3 데이터

[data, time] = loader(gdLogFile);

% total
angletest = 1 : length(time);
veltest = 1 : length(time);
postest = 1 : length(time);

% 배속 시뮬레이션
angletest = 259 : 910;
veltest = 1415 : 2067;
postest = 2571 : 3200;

% A3 시뮬레이션
angletest = 2553 : 3204;
veltest = 3456 : 4108;
postest = 4446 : 5097;


range = angletest;
figure(1)
clf
hold on
grid on
plot(time(range),data.rpy_0(range),'k')
plot(time(range),data.rpdCmd_0(range),'r--')
plot(time(range),islocalmax(data.rpy_0(range))*5+mean(data.rpy_0(range)))
plot(time(range),-islocalmin(data.rpy_0(range))*5+mean(data.rpy_0(range)))
title('angle')
size(find(islocalmax(data.rpy_0(range))),1)

range = veltest;
figure(2)
clf
hold on
grid on
plot(time(range),data.velNed_1(range),'k')
plot(time(range),data.velCmdNav_1(range),'r--')
plot(time(range),islocalmax(data.velNed_1(range))/5+mean(data.velNed_1(range)))
plot(time(range),-islocalmin(data.velNed_1(range))/5+mean(data.velNed_1(range)))
title('vel')
size(find(islocalmax(data.velNed_1(range))),1)

range = postest;
figure(3)
clf
hold on
grid on
plot(time(range),data.posNed_1(range),'k')
plot(time(range),data.posCmdNed_1(range),'r--')
plot(time(range),islocalmax(data.posNed_1(range))/5+mean(data.posNed_1(range)))
plot(time(range),-islocalmin(data.posNed_1(range))/5+mean(data.posNed_1(range)))
% plot(time(range),data.rpy_0(range)/10+mean(data.posNed_1(range)))
% plot(time(range),data.velNed_1(range)+mean(data.posNed_1(range)))
title('pos')
size(find(islocalmax(data.posNed_1(range))),1)

figure(4)
clf
hold on
grid on
title('response time plot')
plot(data.rpdCmd_0,'k')
plot(data.posCmdNed_1,'r')
plot(data.velCmdNav_1,'m')

%% Reference sweep signal

figure(5)
clf
c1 = 4; c2 = 0.01866;
wmin = 0.4*2*pi;
wmax = 10*2*pi;
t = 0:0.02:13;
T = 13;
res = c2*(wmax-wmin)*T/c1;
omega = wmin + (exp(t/T*c1)-1)*c2*(wmax-wmin);
theta = wmin*t + c2*(wmax-wmin)*(T/c1*exp(c1/T*t)-t) - res;
plot(t,sin(theta))
plot(t,omega/2/pi)