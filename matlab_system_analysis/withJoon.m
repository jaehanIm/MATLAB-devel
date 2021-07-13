clear input
% [data,time] = loader('/home/jaehan/log/gdLogCsv/gdLog_210413_150253.csv');

[data,time] = loader('/home/jaehan/log/gdLogCsv/gdLog_210419_120245.csv');


startIdx = 1;
endIdx = size(data,1);

% startIdx = 21032;
% endIdx = 22077;
% 
% startIdx = 19587;
% endIdx = 20237;

figure(1)
clf
grid on
hold on
plot(data.fcMcMode(startIdx:endIdx))
plot(data.posNed_0(startIdx:endIdx)*100)
plot(data.posNed_1(startIdx:endIdx)*100)

figure(2)
clf
grid on
hold on

plot(data.posNed_0(startIdx:endIdx)*100,data.posNed_1(startIdx:endIdx)*100,'*--')
% plot(data.posNed_0(i)*100,data.posNed_1(i)*100,'*--')

axis equal

for i = 1:size(hope2,2)
    input(i,:) = hope2(i).Position;
end

% input(:,1) = data.posNed_0 * 100;
% input(:,2) = data.posNed_1 * 100;

[xc,yc,R,a] = circfit(input(:,1),input(:,2))

figure(3)
clf
grid on
hold on
plot(input(:,1),input(:,2),'b*--')
th = linspace(0,2*pi,100);
plot(xc + R * cos(th), yc + R*sin(th),'*--')
axis equal

