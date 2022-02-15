%% Actual source
gdLog = readtable('/home/jaehan/log/220121_gimbalTest/gdLog_220121_101056.csv');
gdLog = gdLog(8750:end,:); %pitch 
% gdLog = gdLog(6000:end,:); %roll
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);

% select data
Y = gdLog.gimbalRpy_deg_1;
U = gdLog.gimbalRpyCmd_deg_1;

figure(1)
clf
hold on
grid on
plot(gdTimeS,Y)

%% linearize function
L = length(Y);
T = linspace(0,gdTimeS(end),L);

%% tf estimation
fun = @(x) sum(((Y - lsim(feedback(...
    tf([x(2) x(3) x(4)],[1 0])*...
    tf(x(1),[1 0 0]),1),U,T)).^2));

% GA를 쓰고싶다면 아래 두 줄 사용
options = optimoptions('ga','ConstraintTolerance',1e-6);
result = ga(fun,4,[],[],[],[],[0 0 0 0],[inf inf inf inf],[],options);

% fmincon 쓰고싶다면 아래 두 줄 사용
% options = optimoptions('fmincon','ConstraintTolerance',1e-10);
% result = fmincon(fun,[1 1 1 1],[],[],[],[],[0 0 0 0],[inf inf inf inf],[],options);

%% plot result
figure(2)
clf
hold on
grid on
x = result;
lsim(feedback(...
    tf([x(2) x(3) x(4)],[1 0])*...
    tf(x(1),[1 0 0]),1),U,T)
plot(T,Y,'r')
disp(num2str(result))
feedback(...
    tf([x(2) x(3) x(4)],[1 0])*...
    tf(x(1),[1 0 0]),1)

prelog = x;