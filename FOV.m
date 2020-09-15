title('FOV fail rate')
close all
sampleNum = 1e4;
e_yaw = normrnd(0,0.3,[1 sampleNum]);
e_r = normrnd(0,0.5,[1 sampleNum]);
[X,Y] = meshgrid(-3.1:0.05:3.1,-3.1:0.05:3.1);
out = zeros(size(X));
for i = 1:size(X,1)
for j= 1:size(X,2)
px = length(find(abs(e_yaw-X(i,j))<3.1))/length(e_yaw);
py = length(find(abs(e_r-Y(i,j))<3.1))/length(e_r);
out(i,j) = 1-px*py;
end
end
mesh(X,Y,out)
hold on
mesh(X(24:end-24,23:end-23),Y(24:end-24,23:end-23),out(24:end-24,23:end-23),'edgeColor',[0 0.7 0])
title('>99% area')
xlabel('x pos [m]')
ylabel('y pos [m]')
view(0,90)