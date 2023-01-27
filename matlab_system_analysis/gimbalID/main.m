data = readtable('/home/jaehan/log/gimbal_system_ID/test_gimbal_step_50hz.csv');

L = size(data,1);
r2d = 180/pi;
time = 0:0.02:0.02*(L-1);
time = time';
% time = zeros(L,1);
% for i = 1:L
%     time(i) = datenum(data.time{i}(11:end));
% end
% time = (time-time(1));

figure(1)
clf
subplot(2,1,1)
plot(data.CMD_R*r2d)
hold on
plot(data.CUR_R*r2d)
subplot(2,1,2)
plot(data.CMD_P*r2d)
hold on
plot(data.CUR_P*r2d)

figure(2)
clf
subplot(2,1,1)
plot(time,data.CMD_R*r2d)
hold on
plot(time,data.CUR_R*r2d)
subplot(2,1,2)
plot(time,data.CMD_P*r2d)
hold on
plot(time,data.CUR_P*r2d)

%% Parse data
global nD nN
nN = 2;
nD = 3;

% intIdx = [1576, 2041]; % roll test
intIdx = [2198, 2615]; % pitch test

%% TF estimation

% response = data.CUR_R(intIdx(1):intIdx(2));
% cmd = data.CMD_R(intIdx(1):intIdx(2));

response = data.CUR_P(intIdx(1):intIdx(2));
cmd = data.CMD_P(intIdx(1):intIdx(2));

[Num, Den, delay] = estimate_tf(response,cmd);
sys = tf(Num,Den,'ioDelay',delay)

%% System Characteristic
stepinfo(sys)
figure(3)
clf
step(sys)
title('Step response : Pitch')
temp = step(sys);

figure(4)
clf
margin(sys)
title('Bode plot : Pitch')
[gm,pm,wg,wp] = margin(sys)
sse = 1-temp(end)

%% function 
function [Num, Den, delay] = estimate_tf(Response,Cmd)
global nD nN

timeseriesSet = iddata(Response,Cmd,0.02);
sys = tfest(timeseriesSet,nD,nN,nan);
Num = sys.Numerator;
Den = sys.Denominator;
delay = sys.IODelay;

end