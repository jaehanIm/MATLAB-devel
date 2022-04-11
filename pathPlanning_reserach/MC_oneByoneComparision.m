Main
Main_soleACO

disp("================")
degreeConnectivity
                
clusteringTime
interCompleteTime + intraCompleteTime
solveTime
totalScoreHistory

completeTime_soleACO
soleACO_time
soleACO_result

%%
nodeNumX = [154 154 108 476 973 3234 475];
dconnX = [0.05 0.14 0.06 0.02 0.01 0.0061 0.4698];
connTACO = [0.36 0.45 0.15 8.86 69.65 2975.9 9.8581];
solveTACO = [86.03 74.23 42.01 315.69 1071.7 11109 377.4195];
resACO = [1987.8 1745.7 1564.7 3491.2 5229.8 9984.5 3437.7];
clusTME = [0.82 0.77 0.71 1.2 2.83 12.03 1.8482];
connTME = [0.39 0.50 0.22 5.47 43.19 2061.6 18.25];
solveTME = [60 51.62 35.36 109.14 743.23 2336.6 263.3285];
resTME = [1625.6 1504.5 1283.6 2880.4 3958.3 7620.8 2903.7];
figure(1111)
clf
loglog(nodeNumX,solveTACO,'o','LineWidth',2)
hold on
grid on
loglog(nodeNumX,solveTME,'o','LineWidth',2)
% loglog(nodeNumX,solveTACO)
% loglog(nodeNumX,solveTME)
hold off
legend('ACO','Mine')
title('solution time')
xlabel('node number')


figure(2222)
clf
semilogx(nodeNumX,resACO./resTME,'o','LineWidth',2)
title('solution ratio')
xlabel('node number')
grid on

figure(3333)
clf
semilogy(dconnX,resACO./resTME,'o','LineWidth',2)
title('solution ratio')
xlabel('connectivity')
grid on

figure(4444)
clf
semilogx(dconnX,solveTACO,'o','LineWidth',2)
hold on
grid on
semilogx(dconnX,solveTME,'o','LineWidth',2)
title('solution time')
xlabel('connectivity')
grid on