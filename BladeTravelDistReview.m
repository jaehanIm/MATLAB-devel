clear all
close all
for num = 1:8
idx = num;

logList = ["200703_140828","200703_140905","200703_140940"...
    ,"200703_141016","200703_141051","200703_141127","200703_141202","200703_141238"];
% logList = ["200630_161711","200630_161829","200630_161943"];
directory = '/home/jaehan/log/';

path = directory+logList(idx)+'/spdlog_'+logList(idx)+'.txt';
path = convertStringsToChars(path);
fileID = fopen(path);

n = 1;
while ~feof(fileID)
    data{n} = fgetl(fileID);
    n = n+1;
end

remem = [];
for i = 1:size(data,2)
    temp = strfind(data{i},'BTD');
    if isempty(temp)
        remem(i) = nan;
    else
        remem(i) = temp;
    end
end

n=1;
for i = 1:size(data,2)
    if ~isnan(remem(i))
        tipDist(n) = str2num(data{i}(50:end-1));
        n = n+1;
    end
end

result{num} = tipDist;
end

figure()
hold on
for i = 1:8
plot(result{i})
end
