function [data, time] = loader(gdLogFile)

timeFormat = 1; % Time format (1: seconds, 2: minute, 3: date)

gdLog = importfile(gdLogFile);
data = gdLog;
timeValidIdx = max(find(gdLog.rosTime < 1514764800)); % 2018/01/01 == 1514764800
data((1:timeValidIdx),:) = [];
dataSize = size(data.rosTime);

fTimeGlobal = data.rosTime;
fTimeGlobal_Date = datetime(data.rosTime,'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
fTimeLocal = fTimeGlobal(:,1) - fTimeGlobal(1,1);
fTimeLocal_Sec = seconds(fTimeLocal);

time = fTimeLocal_Sec;
timeSec = fTimeLocal_Sec;
hourTime = fTimeGlobal_Date;
dateStart = fTimeGlobal_Date(1);

switch timeFormat % (1: seconds, 2: minute, 3: date)
    case 1
        strTimeXtickFormat = 's';
    case 2
        strTimeXtickFormat = 'mm:ss.SSS';
    case 3
        time = time + dateStart;
        strTimeXtickFormat = 'HH:mm:ss';
    otherwise
end

sTime = min(time);
eTime = max(time);
timeDurTotal = [sTime eTime];
