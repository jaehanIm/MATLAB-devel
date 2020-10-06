% Run after running 'loader.m' or 'vehicle_analysis.m'

z = 5;

% Drone pointing accuracy
figure(1)
hold on
clf
x_data = detrend(data.rpy_2(range),0)*pi/180*50;
y_data = detrend(data.rpy_1(range),0)*pi/180*50;
mu_x = mean(x_data);
mu_y = mean(y_data);
std_x = std(x_data);
std_y = std(y_data);
[X,Y] = meshgrid(linspace(mu_x-std_x*z,mu_x+std_x*z,50),linspace(mu_y-std_y*z,mu_y+std_y*z,50));
out = BiNorm(X,Y,[mu_x;mu_y],cov(x_data,y_data));
hold on
mesh(X,Y,out)
contour3(X,Y,out)
colormap(hot)
view(0,90)

title('Drone Pointing Accuracy (50m distance)','FontSize',14)
xlabel('horizontal heading[m]')
ylabel('vertical heading [m]')
axis equal

% Gimbal pointing accuracy
figure(2)
hold on
clf
x_data = detrend(data.gimbalRPY_2(range),0)*pi/180*50;
y_data = detrend(data.gimbalRPY_0(range),0)*pi/180*50;
mu_x = mean(x_data);
mu_y = mean(y_data);
std_x = std(x_data);
std_y = std(y_data);
[X,Y] = meshgrid(linspace(mu_x-std_x*z,mu_x+std_x*z,50),linspace(mu_y-std_y*z,mu_y+std_y*z,50));
out = BiNorm(X,Y,[mu_x;mu_y],cov(x_data,y_data));
hold on
mesh(X,Y,out)
contour3(X,Y,out)
colormap(hot)
view(0,90)

title('Gimbal Pointing Accuracy (50m distance)','fontsize',14)
xlabel('horizontal pointing [m]')
ylabel('vertical pointing[m]')
axis equal

% Gimbal pointing accuracy in polar coordinate
figure(3)
hold on
clf
x_data = detrend(data.gimbalRPY_2(range),0)*pi/180*50;
y_data = detrend(data.gimbalRPY_0(range),0)*pi/180*50;
mu_x = mean(x_data);
mu_y = mean(y_data);
std_x = std(x_data);
std_y = std(y_data);
[X,Y] = meshgrid(linspace(mu_x-std_x*z,mu_x+std_x*z,50),linspace(mu_y-std_y*z,mu_y+std_y*z,50));
out = BiNorm(X,Y,[mu_x;mu_y],cov(x_data,y_data));
hold on
mesh(X,Y,out)
contour3(X,Y,out)
colormap(hot)
view(0,90)

title('Gimbal Pointing Accuracy (50m distance)','fontsize',14)
xlabel('horizontal pointing [m]')
ylabel('vertical pointing[m]')
axis equal


% Attitude correlation
figure(4)
hold on
clf
x_data = detrend(data.rpy_0(range),0);
y_data = detrend(data.rpy_1(range),0);
mu_x = mean(x_data);
mu_y = mean(y_data);
std_x = std(x_data);
std_y = std(y_data);
[X,Y] = meshgrid(linspace(mu_x-std_x*z,mu_x+std_x*z,50),linspace(mu_y-std_y*z,mu_y+std_y*z,50));
out = BiNorm(X,Y,[mu_x;mu_y],cov(x_data,y_data));
hold on
mesh(X,Y,out)
contour3(X,Y,out)
colormap(hot)
view(0,90)

title('Attitude Correlation','fontsize',14)
xlabel('[deg]')
ylabel('[deg]')
axis equal