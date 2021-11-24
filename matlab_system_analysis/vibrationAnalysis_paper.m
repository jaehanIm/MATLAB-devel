%% Load vibration data
temp = load('wDamper_01.mat');
wdata = temp.Low_g_Acceleration;
wdata = wdata(:,9000:60000);
temp = load('woDamper_01.mat');
wodata = temp.Low_g_Acceleration;
% wodata = wodata(:,10000:50000);

%% select data
data1 = wodata(2,:);
data2 = wdata(2,:);
Fs = length(wdata)/wdata(1,end);

%% SVD analysis - feature extraction
L = length(data1);
N = 10;
Nsq = N^2;
interval = floor(L / Nsq);
STE = zeros(Nsq,1);
for i = 1:Nsq
STE(i) = sum(data1(interval * (i-1) + 1 : interval * i).^2)/interval;
end
A = zeros(N,N);
for i = 1:Nsq
A(i) = STE(i);
end
A = A';
svdResult = svd(A);
figure(1)
hold on
grid on
plot(svdResult,'LineWidth',4)
xlabel('feature #')
ylabel('amplitude')
title('STE-SVD vibration feature extraction')

figure(2)
hold on
grid on
plot(STE,'LineWidth',4)

%% SVD analysis - denoising
test = data1(1:1000);
% test_hankel = hankel(test);
colNum = 50;
testLen = length(test);

test_hankel = zeros(testLen+1-colNum,colNum);
for i = 1:testLen + 1 - colNum
    test_hankel(i,:) = test(i:i+colNum-1);
end

[U,S,V] = svd(test_hankel);
test_regen_full = U * S * V';
rate = 10;
test_regen = U(:,1:rate) * S(1:rate,1:rate) * V(:,1:rate)';

signal_regen_r = test_regen(1,:);
signal_regen_c = test_regen(2:end,end)';
signal_regen = [signal_regen_r,signal_regen_c];

figure(4)
clf
hold on
grid on
plot(test,'LineWidth',1)
% plot(test_regen_full(1,:))
plot(signal_regen(1,:),'LineWidth',1)

[freq1, fftResult1, psdResult1] = data2fftpsd(test,Fs);
[freq2, fftResult2, psdResult2] = data2fftpsd(signal_regen(1,:),Fs);

figure(5)
clf
hold on
grid on
plot(freq1,fftResult1)
plot(freq2,fftResult2)

%% 