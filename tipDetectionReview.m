clear all
% close all

idx = 7;

logList = ["200611_111447","200611_113320","200611_125154","200611_131752"...
    ,"200611_135228","200611_140139","200611_143043"];
directory = '/home/jaehan/Desktop/test flight/';

path = directory+logList(idx)+'/spdlog_'+logList(idx)+'.txt';
path = convertStringsToChars(path);
fileID = fopen(path);

n = 1;
while ~feof(fileID)
    data{n} = fgetl(fileID);
    n = n+1;
end

for i = 1:size(data,2)
    temp = strfind(data{i},'Within');
    if isempty(temp)
        remem(i) = nan;
    else
        remem(i) = temp;
    end
end

n=1;
for i = 1:size(data,2)
    if ~isnan(remem(i))
        tipDist(n) = str2num(data{i}(76:end-1));
        n = n+1;
    end
end

figure()
plot(tipDist)

bl0 = [23.89 22.25 25.13 23.67 24.38];
bl1 = [23.83 25.1 25.13 24.3];
bl2 = [22.97 24.56 24.29];
figure(1)
clf
hold on
bar([mean(bl0) mean(bl1) mean(bl2)],'BarWidth',0.5,'FaceColor','k','FaceAlpha',0.5)
set(gca,'XTickLabel',{'Blade 0','Blade 1','Blade 2'})
errorbar([mean(bl0) mean(bl1) mean(bl2)],[std(bl0) std(bl1) std(bl2)],'LineWidth',2,'Color','r')
grid on
plot([0 4],[23 23],'k--','LineWidth',1)
title('Blade Length','fontsize',14)
xlabel('Blade number','fontsize',14)
ylabel('Blade length [m]','fontsize',14)