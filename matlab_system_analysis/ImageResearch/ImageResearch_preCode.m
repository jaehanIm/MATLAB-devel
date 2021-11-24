addpath('..');
% sampleImage = imread('/home/jaehan/Desktop/MATLAB devel/matlab_system_analysis/ImageResearch/BBB_5350.JPG');
% sampleImage = imread('/home/jaehan/Desktop/MATLAB devel/matlab_system_analysis/ImageResearch/BBB_5537.JPG');
sampleImage = imread('/home/jaehan/log/50_set/1104_12s/BBB_5513.JPG'); %흐림
% sampleImage = imread('/home/jaehan/log/50_set/1104_12s/BBB_5514.JPG'); % 선명

sampleImage = im2gray(sampleImage);
sampleImage_d = double(sampleImage);
sampleImage_delX_d = abs(diff(sampleImage_d,1,1));
sampleImage_delY_d = abs(diff(sampleImage_d,1,2));

figure(1)
imshow(sampleImage)

% sampleImage_del_d = sqrt(sampleImage_delX_d(:,1:end-1).^2 + sampleImage_delY_d(2:end,:).^2);
sampleImage_del_d = sampleImage_delX_d(:,1:end-1) .* sampleImage_delY_d(2:end,:);

%% Edge detection

maxDel = max(max(sampleImage_del_d));
normCoeff = 255 / maxDel;

sampleImage_del_d = sampleImage_del_d*normCoeff;
sampleImage_del = uint8(sampleImage_del_d);
figure(2)
imshow(sampleImage_del)
max(sampleImage_delX_d,[],'all')
max(sampleImage_delY_d,[],'all')

% [freq,fftR] = data2fftpsd(im2double(sum(sampleImage_del)),1);
% [freq,fftR] = data2fftpsd(im2double(sampleImage_del(:)),1);

figure(3)
g = histogram(sampleImage_del_d);
g.BinWidth = 2;

figure(4)
thres = 1;
sampleImage_del_thres = sampleImage_del;
sampleImage_del_thres(sampleImage_del_thres<=thres) = 0;
imagesc(sampleImage_del_thres)


%% SVD denoising

% [U,S,V] = svd(sampleImage_del_d);
% 
% % Reconstruction
% 
% regenLength = 200;
% sampleImage_regen_d = U(:,1:regenLength) * S(1:regenLength,1:regenLength) * V(:,1:regenLength)';
% sampleImage_regen = uint8(sampleImage_regen_d);
% 
% figure(4)
% imshow(sampleImage_regen)
% 
% figure(5)
% plot(freq,fftR/length(freq))

%% Partial Analysis
% sharp = sampleImage_del_d(4965,5300: 6000);
% smooth = sampleImage_del_d(4700:5400,5655);
% sharp = sampleImage_del_d(4965:5100,5300:6000);
% smooth = sampleImage_del_d(4700:5400,5655);
% sharp = sampleImage_del_d(4965,:);
% smooth = sampleImage_del_d(:,5655);

% sharp = sharp.^2;
% smooth = smooth.^2;

% figure()
% plot(sharp)
% hold on
% plot(smooth)
% 
% kernelSize = 10;
% xnum = floor(size(sampleImage,1)/kernelSize);
% ynum = floor(size(sampleImage,2)/kernelSize);
% result = zeros(xnum,ynum);
% for i = 1:xnum
%     for j = 1:ynum
%         section = sampleImage_del_d(kernelSize*(i-1)+1:kernelSize*i,kernelSize*(j-1)+1:kernelSize*j);
%         result(i,j) = kurtosis(section,1,'all');
%     end
% end
% 
% result = result/max(result,[],'all')*255;
% result = uint8(result);

% [U,S,V] = svd(hankel(sharp));
% filter = 200;
% sharp_regen = U(:,1:filter) * S(1:filter,1:filter) * V(:,1:filter)';
% sharp_regen = sharp_regen(1,:);
% % [U,S,V] = svd(hankel(smooth));
% 
% figure()
% plot(sharp)
% hold on
% plot(sharp_regen)