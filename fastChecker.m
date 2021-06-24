function fastChecker(logFile)

[data,time] = loader(logFile);

figure(1)
clf
grid on
hold on
plot(time,data.fcMcMode)

end