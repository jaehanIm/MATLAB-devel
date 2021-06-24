title('FOV fail rate')
close all
sampleNum = 1e4;
e_yaw = normrnd(0,std_x,[1 sampleNum]);
e_r = normrnd(0,std_y,[1 sampleNum]);
[X,Y] = meshgrid(-0.92:0.01:0.92,-0.69:0.01:0.69);
out = zeros(size(X));
for i = 1:size(X,1)
for j= 1:size(X,2)
px = length(find(abs(e_yaw-X(i,j))<0.92))/length(e_yaw);
py = length(find(abs(e_r-Y(i,j))<0.69))/length(e_r);
out(i,j) = 1-px*py;
end
end
figure(2)
clf
mesh(X,Y,out-0.1)
hold on
contour3(X,Y,out,'ShowText','on')
title('>99% area')
xlabel('x pos [m]')
ylabel('y pos [m]')
view(-45,60)