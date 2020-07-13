%% data log explorer
total_path = "/home/jaehan/Desktop/test flight/WT_yaw/";
flight_list = dir(total_path);
path_name_list = [];
for i = 3:length(flight_list)
    path_name_list{i-2} = strcat(total_path,flight_list(i).name,'/');
end
result = [];
record_size = [];

for k = 1:length(path_name_list)
    %% data loading
%     path_name = "/home/jaehan/Desktop/test flight/WT_yaw/log_0109_5/";
    path_name = path_name_list{k};
    list = dir(path_name);
    n = length(list);
    file_name = [];

    for i = 3:n
        file_name{i-2} = strcat(path_name,list(i).name);
    end

    %% data extraction
    nf = length(file_name);
    buffer_time = 10; %sec
    buffer = buffer_time * 50;
    % init_yaw = zeros(nf,1);
    c=1;
    init_yaw = [];
    init_yaw_var = [];
    local_result = [];
    
    for i = 1:nf    
        if size(importfile(file_name{i}),1) ~= 0
            [time, out] = logReader(convertStringsToChars(file_name{i}));
            init_time = find(out(:,1) == 1,1);
            if ~isempty(init_time) && init_time > buffer
                init_yaw(c) = mean(out(init_time-buffer:init_time,4));
                init_yaw_var(c) = var(out(init_time-buffer:init_time,4));
                c=c+1;
            end
        end
        disp(['Complete by ',num2str(i),'/',num2str(nf)]);
    end

    local_result = [init_yaw',init_yaw_var'];
    record_size = [record_size,length(local_result)];
    result = [result;local_result];
end

%% log file record
plot(result(:,1),'b:')
hold on
plot(result(:,1),'bo')
% plot(11,result(11,1),'rx','MarkerSize',10,'LineWidth',3)
plot([10.5 10.5],[0 100],'r--')
plot([23.5 23.5],[0 100],'r--')
plot([32.5 32.5],[0 100],'r--')
plot([43.5 43.5],[0 100],'r--')
xlim([0 54])
xlabel('data number')
xlabel('data number','fontsize',14)
ylabel('angle [deg]','fontsize',14)

%% variance 
% start = [1 12 24 33 44];
% ending = [10 23 32 43 53];
start = [5 13 25 34 45];
ending = [10 19 33 41 53];

daily_mean = [];
daily_var = [];
daily_max = [];
daily_min = [];
c = 1;
for k = 1:length(start)
    i = start(k); j = ending(k);
    daily_mean(c) = mean(result(i:j,1));
    daily_var(c) = var(result(i:j,1));
    daily_max(c) = max(result(i:j,1));
    daily_min(c) = min(result(i:j,1));
    c=c+1;
end

figure()
hold on
plot(daily_max-daily_mean,'r*')
plot(daily_min-daily_mean,'b*')
plot(daily_var,'r:')
plot(-daily_var,'r:')
plot([1 5],[0 0],'k-')
grid on

%% season
summer = []; winter = []; angle_data = [];
summer = vertcat(result(start(2):ending(2),1)-daily_mean(2),result(start(3):ending(3),1)-daily_mean(3));
winter = vertcat(result(start(1):ending(1),1)-daily_mean(1),result(start(4):ending(4),1)-daily_mean(4),result(start(5):ending(5),1)-daily_mean(5));

% season = [];
for i = 1:length(summer)
    season(i,:) = char('summer');
end
for i = length(summer)+1:length(summer)+length(winter)
    season(i,:) = char('winter');
end

angle_data = vertcat(summer,winter);

figure()
boxplot(angle_data(:,1),season);

figure()
histogram(summer,'facealpha',.5,'Normalization','probability')
hold on
grid on
histogram(winter,'FaceAlpha',.5,'Normalization','probability')

