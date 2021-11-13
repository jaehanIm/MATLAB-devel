image_real = imread('/home/jaehan/Desktop/MATLAB devel/Photo_object_detection/TTT_0818.JPG'); % clear
% image_real = imread('/home/jaehan/Desktop/MATLAB devel/Photo_object_detection/TTT_1265.JPG'); % scattered
% image_real = imread('/home/jaehan/Desktop/MATLAB devel/Photo_object_detection/TTT_0456.JPG'); % bad weather
% image_real = imread('/home/jaehan/Desktop/MATLAB devel/Photo_object_detection/TTT_0457.JPG'); % dark 

figure(1)
histogram(image_real(:,:,1),'Normalization',"count")
grid on
figure(2)
histogram(image_real(:,:,2),'Normalization',"count")
grid on
figure(3)
histogram(image_real(:,:,3),'Normalization',"count")
grid on

image_sky = image_real(1:end,100:200,:);
figure(8)
imshow(image_sky)
figure(9)
histogram(image_sky(:,:,1),'Normalization',"count")
grid on
figure(10)
histogram(image_sky(:,:,2),'Normalization',"count")
grid on
figure(11)
histogram(image_sky(:,:,3),'Normalization',"count")
grid on

image_blade = image_real(1:end,4300:4400,:);
figure(4)
imshow(image_blade)
figure(5)
histogram(image_blade(:,:,1),'Normalization',"count")
grid on
figure(6)
histogram(image_blade(:,:,2),'Normalization',"count")
grid on
figure(7)
histogram(image_blade(:,:,3),'Normalization',"count")
grid on

image_gray = rgb2gray(image_real);
figure(12)
imshow(image_gray)
figure(13)
histogram(image_gray,'Normalization',"count")
grid on

figure(14)
imshow(image_real)
figure(15)
temp = diff(image_gray) == 0;
imshow(temp)
figure(16)
imshow(diff(image_real).^10)