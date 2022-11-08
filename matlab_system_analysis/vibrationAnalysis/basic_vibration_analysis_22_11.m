%% Load and parse data
% data = readtable('/home/jaehan/log/221108_153535/aSensorImu_221108_153535.csv');
data = readtable('/home/jaehan/log/221108_154600/aSensorImu_221108_154600.csv');

for i = 1:size(data,1)
    if data.rosTime(i) > 100000
        thres = i;
        initTime = data.rosTime(i);
        break;
    end
end

%% raw data plot
figure(1)
clf
plot(data.rosTime(thres:end)-initTime,data.acc_mpss_0(thres:end))
hold on
plot(data.rosTime(thres:end)-initTime,data.acc_mpss_1(thres:end))
plot(data.rosTime(thres:end)-initTime,data.acc_mpss_2(thres:end))
grid on
figure(2)
clf
plot(data.rosTime(thres:end)-initTime,data.gyro_dps_0(thres:end))
hold on
plot(data.rosTime(thres:end)-initTime,data.gyro_dps_1(thres:end))
plot(data.rosTime(thres:end)-initTime,data.gyro_dps_2(thres:end))
grid on