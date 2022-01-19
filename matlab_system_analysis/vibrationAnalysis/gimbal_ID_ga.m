%% example source
k = 1; D = 5; P = 10; I = 0.2;
sys = tf([k*D k*P k*I],[1 k*D k*P k*I]);
[Y,T] = step(sys);
L = length(Y);
Y = Y + (randn(L,1)-0.5)/10000;

%% linearize function
L = length(Y);
T = linspace(0,T(end),L);

%% tf estimation
fun = @(x) sum((Y - lsim(tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)]),ones(1,L),T)).^4);
% fun = @(x) sum((Y - lsim(tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)]),U,T)).^2);
options = optimoptions('ga','ConstraintTolerance',1e-10);
result = ga(fun,4,[],[],[],[],[0 0 0 0],[],[],options);
% result = fmincon(fun,[0 0 0 0],[],[],[],[],[0 0 0 0],[],[],options);

%% plotting
figure(1)
clf
hold on
grid on
x = result;
lsim(tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)]),ones(1,L),T)
plot(T,Y,'r')

disp('original')
textarr = [num2str(k),' ', num2str(D),' ', num2str(P),' ', num2str(I)];
disp(textarr)
disp(num2str(result))

tf([x(1)*x(2) x(1)*x(3) x(1)*x(4)],[1 x(1)*x(2) x(1)*x(3) x(1)*x(4)])