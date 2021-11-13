sampleImage = imread('/home/jaehan/Desktop/MATLAB devel/matlab_system_analysis/ImageResearch/BBB_5350.JPG');
% sampleImage = imread('/home/jaehan/Desktop/MATLAB devel/matlab_system_analysis/ImageResearch/BBB_5537.JPG');
sampleImage = im2gray(sampleImage);
sampleImage_delX = diff(sampleImage,1,1);
sampleImage_delY = diff(sampleImage,1,2);

sampleImage_del = sampleImage_delX(:,1:end-1) + sampleImage_delY(1:end-1,:);

maxDel = max(max(sampleImage_del));
normCoeff = 255 / maxDel;

sampleImage_del = sampleImage_del * normCoeff;

figure(1)
imshow(sampleImage_del)


%% SVD denoising
[U,S,V] = svd(im2double(sampleImage_del));

%% Reconstruction

regenLength = 100;
sampleImage_regen = U(:,1:regenLength) * S(1:regenLength,1:regenLength) * V(:,1:regenLength)';

figure(2)
imshow(sampleImage_regen)