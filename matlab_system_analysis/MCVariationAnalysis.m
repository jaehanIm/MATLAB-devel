%% Load Mission
[data,time] = loader('/home/jaehan/log/gdLogCsv/gdLog_210309_155929.csv');   % put gdLog address

%% Parameter setting
confidence = 0.9;       % set the value of confidence range to be highlighted [0~1]
density = 500;          % mvncdf calculation cell resolution / adjust this value to control computation time
range = 14508 : 15602;  % specify the range of data index.

%% Internal parameters // DISREGARD
level = 1-confidence;
zRange = 5;
cumulated = 0;

index = [];
peripheral = [];

%% Data selection // Change if required

% Pointing Accuracy - example set
dataset = [posXyz_0(range),posXyz_1(range),posXyz_2(range),data.gimbalRPY_1(range)*d2r,data.gimbalRPY_2(range)*d2r];
num = 1e5;
covData = cov(dataset);
randomData = mvnrnd(zeros(size(dataset,2),1),covData,num);
l_err = 25 + randomData(:,1);
p_err_x = l_err.*tan(randomData(:,5))+randomData(:,2);
p_err_y = l_err.*tan(randomData(:,4))+randomData(:,3);
x_data = p_err_x;
y_data = p_err_y;

%% Main
covar = cov([x_data,y_data]);

x_data = detrend(x_data,0);
y_data = detrend(y_data,0);

mu_x = mean(x_data);
mu_y = mean(y_data);
std_x = std(x_data);
std_y = std(y_data);
maxpdf = mvnpdf([mu_x mu_y],[mu_x mu_y],covar);
thres = level * maxpdf;

xGrid = linspace(mu_x-zRange*std_x,mu_x+zRange*std_x,density);
yGrid = linspace(mu_y-zRange*std_y,mu_y+zRange*std_y,density);
xGrid(end) = [];
yGrid(end) = [];
xGridSize = xGrid(2)-xGrid(1);
yGridSize = yGrid(2)-yGrid(1);
gridArea = xGridSize * yGridSize;

[X,Y] = meshgrid(linspace(mu_x-std_x*zRange,mu_x+std_x*zRange,50),linspace(mu_y-std_y*zRange,mu_y+std_y*zRange,50));
out = BiNorm(X,Y,[mu_x;mu_y],cov(x_data,y_data));

for i = xGrid
    for j = yGrid
        if mvnpdf([i j],[mu_x mu_y],covar) >= thres
            cumulated = cumulated + mvncdf([i j],[i+xGridSize j+yGridSize],[mu_x mu_y],covar);
            index = vertcat(index,[i,j]);
        end
    end
end

cumulated % for verification. This value must be equal to 'confidence' param.

for i = 1:length(index)
    if boxed(i,index,xGridSize,yGridSize)
        peripheral = vertcat(peripheral,index(i,:));
    end
end

figure(1)
clf
hold on
plot(x_data,y_data,'.')
contour(X,Y,out/maxpdf,[0.5 0.2 0.1 0.05 0.01],'ShowText','on')
grid on
plot(peripheral(:,1),peripheral(:,2),'.r','MarkerSize',10)
axis equal
title('confidence area')

confidence_area = length(index) * gridArea;

A = [num2str(confidence*100),'% confidence area is ',num2str(confidence_area),' [deg^2]'];
disp(A);

function flag = boxed(n,index,xGridSize,yGridSize)
x = index(n,1); y = index(n,2);
right = 0; left = 0; up = 0; down = 0;
flag = 0;

if isempty(find(abs(index(:,1)-(x+xGridSize)) < xGridSize/2 & abs(index(:,2)-y) < yGridSize/2,1))
    right = 1;
end
if isempty(find(abs(index(:,1)-(x-xGridSize)) < xGridSize/2 & abs(index(:,2)-y) < yGridSize/2,1))
    left = 1;
end
if isempty(find(abs(index(:,1)-x) < xGridSize/2 & abs(index(:,2)-(y+yGridSize)) < yGridSize/2,1))
    up = 1;
end
if isempty(find(abs(index(:,1)-x) < xGridSize/2 & abs(index(:,2)-(y-yGridSize)) < yGridSize/2,1))
    down = 1;
end

if up || down || right || left
    flag = 1;
end

end