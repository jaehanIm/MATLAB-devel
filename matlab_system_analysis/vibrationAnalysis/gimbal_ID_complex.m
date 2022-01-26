%% Actual source
gdLog = readtable('/home/jaehan/log/220121_gimbalTest/gdLog_220121_101056.csv');
% gdLog = gdLog(8750:end,:);%pitch
gdLog = gdLog(6000:end,:);%roll
gdTimeS = gdLog.rosTime - gdLog.rosTime(1);

% select data
Y = gdLog.gimbalRpy_deg_0;
U = gdLog.gimbalRpyCmd_deg_0;

%% example source
% k = 1; D = 5000; P = 100; I = 2; zeta = 16; omega = 50;
% sys = tf([k*D k*P k*I],[1 2*zeta*omega omega^2 k*D k*P k*I]);
% sys = feedback(tf([k*D k*P k*I],[1 2*zeta*omega omega^2 0 0 0]),1);
% [Y,T] = step(sys);
% L = length(Y);
% Y = Y + (randn(L,1)-0.5)/10000;
% U = ones(L,1);
% figure(1)
% step(sys)


%% linearize function
L = length(Y);
T = linspace(0,gdTimeS(end),L);

%% tf estimation
% x - k D P I zeta omega ioDelay

% fun = @(x) sum(((Y - lsim(...
%     tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 2*x(5)*x(6) x(6)^2 x(1)*x(2) x(1)*x(3) x(1)*x(4)])...
%     ,U,T)).^2));

fun = @(x) sum(((Y - lsim(feedback(...
    tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 2*x(5)*x(6) x(6)^2 0 0 0]),1)...
    ,U,T)).^2));

% GA를 쓰고싶다면 아래 두 줄 사용
% options = optimoptions('ga','ConstraintTolerance',1e-3);
% result = ga(fun,6,[],[],[],[],[0 0 0 0 0 0],[],[],options);

% fmincon 쓰고싶다면 아래 두 줄 사용
options = optimoptions('fmincon','ConstraintTolerance',1e-8);
[result,score] = fmincon(fun,[prelog(1) prelog(2) prelog(3) prelog(4) 1 100 0],[],[],[],[],[],[],[],options);

%% plot result
figure(2)
clf
hold on
grid on
x = result
openTf = tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 2*x(5)*x(6) x(6)^2 0 0 0],'ioDelay',x(7));
closedTf = feedback(openTf,1);
lsim(closedTf,U,T);
plot(T,Y,'r')
% disp(num2str(result))
% [resultNum,resultDen] = ss2tf(closedTf.A,closedTf.B,closedTf.C,closedTf.D);
% resultingSys = tf(resultNum,resultDen,'ioDelay',closedTf.InternalDelay);
resultingSys = closedTf
openTf
score

figure(3)
bode(resultingSys)

figure(33)
bode(openTf)

figure(4)
step(resultingSys)