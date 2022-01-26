%% Actual source
gdLog = readtable('/home/jaehan/log/220121_gimbalTest/gdLog_220121_101056.csv');
% gdLog = gdLog(8750:end,:); %pitch 
gdLog = gdLog(6000:end,:); %roll
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);

% select data
Y = gdLog.gimbalRpy_deg_0;
U = gdLog.gimbalRpyCmd_deg_0;

figure(1)
clf
hold on
grid on
plot(gdTimeS,Y)

%% linearize function
L = length(Y);
T = linspace(0,gdTimeS(end),L);

%% tf estimation
fun = @(x) sum(((Y - lsim(tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)]),U,T)).^2));

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
lsim(tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)]),U,T)
plot(T,Y,'r')
disp(num2str(result))
tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)])

prelog = x;