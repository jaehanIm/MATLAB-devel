%% example source
k = 1; D = 5; P = -0.1; I = 0;
sys = tf([k*D k*P k*I],[1 k*D k*P k*I]);
[Y,T] = step(sys);
L = length(Y);
Y = Y + (randn(L,1)-0.5)/100;

figure(2)
clf
step(sys)

%% equi-interval
L = length(Y);
T = linspace(0,T(end),L);

%% tf estimation
fun = @(x) sum((Y - lsim(tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)]),ones(1,L),T)).^2);

options = optimoptions('ga','ConstraintTolerance',1e-4);
result = ga(fun,4,[],[],[],[],[0 0 0 0],[],[],options);

% options = optimoptions('fmincon','ConstraintTolerance',1e-10);
% result = fmincon(fun,[1 3 10 1],[],[],[],[],[0 -inf -inf -inf],[inf inf inf inf],[],options);

%% plot result
figure(1)
clf
hold on
grid on
x = result;
% estimated tf response
lsim(tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)]),ones(1,L),T)
% actual response
plot(T,Y,'r')

disp('original')
textarr = [num2str(k),' ', num2str(D),' ', num2str(P),' ', num2str(I)];
disp(textarr)
disp(num2str(result))

tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)])