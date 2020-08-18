d2r = pi/180; r2d = 1/d2r;
%% Parameter Setting

battTime = 10; % effective minute for mission
battmarg = 10; % battery life safety margin [%]
int_wtSize = 35; % blade length of wt of interest
int_availft = 500; % available flight time
coc_thres = 0.025; % CoC threshold
std_e_yaw = 0.3116 * d2r; % rad std.
std_e_r = 0.9206; % m std.
%% Data Loading

time6 = load('time6mission.mat');
time4 = load('time4mission.mat');
time6_1 = time6.time6_1;
time6_2 = time6.time6_2;
time4_1 = time4.data;

tot6_1 = time6_1(:,1);
tot6_2 = time6_2(:,1);
tot4_1 = time4_1(:,1);

bat_safeline = (battTime*(1-battmarg/100));

%% Mission Time Composition

data(:,:,1) = time4_1(:,2:end);
data(:,:,2) = time6_1(:,2:end);
data(:,:,3) = time6_2(:,2:end);
data = permute(data,[1 3 2]);
plotBarStackGroups(data,{'23m','35m','49m','70m'})
clc
hold on
plot([1,2,3,4],time6_1(:,1),'k')
plot([1,2,3,4]-0.217,time4_1(:,1),'k')
plot([1,2,3,4]+0.217,time6_2(:,1),'k')
plot([0,5],ones(1,2)*battTime*60,'k-');                       % battTime limit
plot([0,5],ones(1,2)*(battTime*(1-battmarg/100))*60,'k:');   
text(1,battTime*60-20,'Max flight time','FontWeight','bold')
grid on;
title('Job Time Composition','FontSize',14)
xlabel('WT size')
ylabel('Mission time [s]')

disp('===== MTC analysis =====');
temp1 = time6_1(:,2:end)./sum(time6_1(:,2:end),2); temp1 = mean(temp1) * 100; comp1 = temp1;
T1 = ['MTC for ','mission 6_1 : ',num2str(temp1(1)),'% ',num2str(temp1(2)),'% ',num2str(temp1(3)),'% '];
temp2 = time6_2(:,2:end)./sum(time6_2(:,2:end),2); temp2 = mean(temp2) * 100; comp2 = temp2;
T2 = ['MTC for ','mission 6_2 : ',num2str(temp2(1)),'% ',num2str(temp2(2)),'% ',num2str(temp2(3)),'% '];
temp3 = time4_1(:,2:end)./sum(time4_1(:,2:end),2); temp3 = mean(temp3) * 100; comp3 = temp3;
T3 = ['MTC for ','mission 4_1 : ',num2str(temp3(1)),'% ',num2str(temp3(2)),'% ',num2str(temp3(3)),'% '];
disp(T1);disp(T2);disp(T3);
disp(' ');

temp1 = time6_1(:,2:end)./sum(time6_1(:,2:end),2); temp1 = temp1 * 100;
temp2 = time6_2(:,2:end)./sum(time6_2(:,2:end),2); temp2 = temp2 * 100;
temp3 = time4_1(:,2:end)./sum(time4_1(:,2:end),2); temp3 = temp3 * 100;
%% Service Coverage Curve

wtSize = [0.672, 1, sqrt(2), 2] * 35.7;
xrange = wtSize(1)-20:wtSize(end)+20;

f_1 = polyfit(wtSize,tot6_1,2);
f_2 = polyfit(wtSize,tot6_2,2);
f_3 = polyfit(wtSize,tot4_1,2);
y_1 = polyval(f_1,xrange);
y_2 = polyval(f_2,xrange);
y_3 = polyval(f_3,xrange);

figure(2)
clf
grid on
hold on
xlabel('WT size [m]')
ylabel('Total mission time [s]')
title('Service Coverage Curve','fontsize',14)
plot(wtSize,tot6_1,'bo')                                                      % tot6_1
plot(xrange,y_1,'b--')
plot(wtSize,tot6_2,'ro')                                                      % tot6_2
plot(xrange,y_2,'r--')
plot(wtSize,tot4_1,'ko')
plot(xrange,y_3,'k--')
plot(xrange,ones(1,length(xrange))*battTime*60,'k-');                       % battTime limit
bat_safeline = (battTime*(1-battmarg/100));
plot(xrange,ones(1,length(xrange))*bat_safeline*60,'k:');             % battTime safety margin
% plot([int_wtSize int_wtSize],[0 600],'m--')
text(min(xrange)+1,battTime*60+20,'Max flight time','FontWeight','bold')
% ylim([0 max(y_1)+100])
xlim([xrange(1) xrange(end)])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
disp('===== SCC analysis =====');
T = ['Service target wt size : ',num2str(int_wtSize),'m'];
required_ft = polyval(f_1,int_wtSize);
T2 = ['Estimated service time : ',num2str(required_ft),'s'];
T3 = ['Service margin : ',num2str(bat_safeline*60-required_ft),'s'];
avail_wtSize = quadinv(f_1,int_availft);
T4 = ['Maximum service capable wt size : ',num2str(avail_wtSize(1)),'m'];


disp(T); disp(T2); disp(T3); disp(T4);
disp(' ');

%% Image Quality Index
disp('===== Image Quality analysis =====')
sampleNum = 1e5;
disp(['Monte Carlo simulation sample size : ',num2str(sampleNum)])

e_yaw = normrnd(0,std_e_yaw,[1 sampleNum]);
e_r = normrnd(0,std_e_r,[1 sampleNum]);
S2 = (7+e_r)./cos(e_yaw);

S1 = 7; %target focus distance [m]
N  = 3.5; % f#
f  = 43; % focal length [mm]
S2 = S2 * 1000;
S1 = S1 * 1000;
coc = abs(S2-S1)./S2.*f^2./(N.*(S1-f));

figure(3)
clf
grid on; hold on;
h = histogram(coc,'Normalization','probability','FaceColor',[0.2 0.4470 0.7410]*0.5);
plot([0.024 0.024],[0 max(h.Values)],'k--')
portion = length(find(h.Data<coc_thres))/length(h.Data)*100;
T1 = ['Probability of satisfying IQ criteria : ',num2str(portion),' %'];
title('CoC distribution (Hovering)', 'FontSize', 14)
xlabel('CoC [mm]','FontSize',14)
ylabel('%','fontsize',14)
disp(T1);
disp(' ');

%% Sensitivity analysis

disp('===== Sensitivity analysis =====');
slope = (polyval(f_1,int_wtSize+0.1) - polyval(f_1,int_wtSize-0.1))/0.2;
MTC6_1 = temp1./100;
compSlope = MTC6_1.*(1/slope*60);
T1 = ['1. Performance slope : ',num2str(slope),' s/m  or  ',num2str(1/slope*60),' m/min'];
T2 = ['2. MTC performance sensitivity [m/min]']; T4 = num2str(compSlope);
T3 = ['  std       ','  cam       ','  wpt     '];
disp(T1); disp(T2); disp(T3); disp(T4); 
disp(' ');

temp = load('cocPassRatio.mat');
a = temp.a; b = temp.b; portion = temp.portion;
[X,Y] = meshgrid(a,b);
out = quadResponseSurface(a,b,portion);
R = out(1).*X.^2 + out(2).*Y.^2 + out(3).*X.*Y + out(4).*X + out(5).*Y + out(6);
R_exact = out(1)*std_e_yaw^2 + out(2)*std_e_r^2 + out(3)*std_e_yaw*std_e_r + out(4)*std_e_yaw + out(5)*std_e_r + out(6);
figure(4)
clf
mesh(b,a*r2d,R)
xlabel('Distance hold std.dev.');
ylabel('Yaw angle hold std.dev.');
title('CoC Pass ratio to error','FontSize',14)

disp('3. Image quality sensitivity')
T1 = ['CoC sensitivity of std.dev.(e_r) : ',num2str(2*out(1)*std_e_yaw+out(3)*std_e_r+out(4)),'mm/m'];
T2 = ['CoC sensitivity of std.dev.(e_yaw) : ',num2str(2*out(2)*std_e_r+out(3)*std_e_yaw+out(5)),'mm/rad'];
% T1 = ['CoC sensitivity of std.dev.(e_r) : ',num2str(2*out(1)*std_e_r+out(3)*std_e_yaw+out(4)),'mm/m'];
% T2 = ['CoC sensitivity of std.dev.(e_yaw) : ',num2str(2*out(2)*std_e_yaw+out(3)*std_e_r+out(5)),'mm/rad'];
disp(T1);
disp(T2);