%% Variable setting

S1 = 7; %target focus distance [m]
N  = 3.5; % f#
f  = 43; % focal length [mm]

%% Data loading
% gdLog = '/home/jaehan/Desktop/test flight/gdLogCsv/gdLog_200414_102453.csv';
gdLog = '/home/jaehan/Desktop/test flight/gdLogCsv/gdLog_200414_113807.csv';
% gdLog = '/home/jaehan/Desktop/test flight/gdLog_200520_153847.csv';

[time, out, acc, gimbal] = logReader(gdLog);

%% Main
S2 = 5:0.01:12;
S2 = S2 * 1000;
S1 = S1 * 1000;
coc = abs(S2-S1)./S2.*f^2./(N.*(S1-f));

figure(1)
hold on
plot(S2/1000,coc)
plot([7 7],[0 0.05],'k--')
plot([5 12],[0.025 0.025],'k--')
grid on
xlabel('Object distance [m]','fontsize',14)
ylabel('CoC (Circle of Confusion) [mm]','fontsize',14)
title('CoC & Error distribution','fontsize',14)

%% Statistics
% sample1 = [52810 53180]; sample2 = [52250 53180];
% sample1 = [15430 18150]; sample2 = [14700 18150];
sample1 = [15430 16000]; sample2 = [14700 15500];
yaw = out(:,5);
north = out(:,6);
east = out(:,7);
height = out(:,8);
ey_=std(yaw(sample1(1):sample1(2)));   % yaw angle std during TE mission
en_=std(north(sample2(1):sample2(2))); % north position std during TE mission
ee_=std(east(sample2(1):sample2(2)));  % east position std during TE mission
eh_=std(height(sample2(1):sample2(2)));% algitude std during TE mission

%% Image quality function (CoC)
num = 1e6;
theta = 45*pi/180;
etheta = normrnd(0,ey_,[1 num])*pi/180;
en = normrnd(0,en_,[1 num]);
ee = normrnd(0,ee_,[1 num]);
eh = normrnd(0,eh_,[1 num]);
temp = ((7+en).*tan(theta+etheta)).^2+(7+en).^2;
temp = sqrt(temp);
temp = sqrt(temp.^2+eh.^2);
% for i = 1:length(theta)
%     temp = ((7+en).*tan(theta(i)+etheta)).^2+(7+en).^2;
%     temp = sqrt(fd);
%     fd(i) = std(temp)
% end
hold on
histogram(temp,'normalization','probability')
xlim([5 12])
ylim([0 0.05])

temp = ((7+en.*0).*tan(theta+etheta)).^2+(7+en.*0).^2;
temp = sqrt(temp);
temp = sqrt(temp.^2+eh.^2);
histogram(temp,'normalization','probability')

temp = 

%% Image quality function (Film angle)
temp2 = etheta + ee./(7+en); %angle
temp3 = temp2 .* (7+en); %distance

%% Image quality history (CoC)
ey_=(yaw(sample2(1):sample2(2)));   
ey_ = ey_-ey_(1);
en_=(north(sample2(1):sample2(2)));
en_ = en_-mean(en_);
ee_=(east(sample2(1):sample2(2)));
ee_ = ee_-mean(ee_);
eh_=(height(sample2(1):sample2(2)));
eh_ = eh_-mean(eh_);

temp = ((7+en_).^2.*tan(ey_*pi/180).^2+(7+en_).^2);
temp = sqrt(temp);
temp_angle = ((7).^2.*tan(ey_*pi/180).^2+(7).^2);
temp_angle = sqrt(temp_angle);
temp = sqrt(temp.^2+(eh_).^2);

figure()
plot(temp)
hold on
plot(temp_angle);
title('Object distance plot','fontsize',14)
xlabel('data #','fontsize',14)
ylabel('dist [m]','fontsize',14)
plot([0 1000],[7 7],'k')
grid on

figure()
S2 = temp*1000;
coc = abs(S2-S1)./S2.*f^2./(N.*(S1-f));
S2 = temp_angle*1000;
coc_angle = abs(S2-S1)./S2.*f^2./(N.*(S1-f));
plot(coc)
grid on
hold on
plot(coc_angle)
plot([0 1000],[0.025 0.025],'k')
ylim([0 0.030])
title('CoC plot','fontsize',14)
xlabel('data #','fontsize',14)
ylabel('CoC [mm]','fontsize',14)