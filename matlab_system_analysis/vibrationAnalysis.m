%% Load vibration data
temp = load('wDamper_01.mat');
wdata = temp.Low_g_Acceleration;
wdata = wdata(:,9000:60000);
temp = load('woDamper_01.mat');
wodata = temp.Low_g_Acceleration;
% wodata = wodata(:,10000:50000);

for i = 1:3
    wdata(i+1,:) = detrend(wdata(i+1,:),1);
    wodata(i+1,:) = detrend(wodata(i+1,:),1);
end

Fs = length(wdata)/wdata(1,end);

figure(1)
title('FFT')
clf
[freq, fftResult, psdResult] = data2fftpsd(wdata(2,:),Fs);
plot(freq,fftResult);
hold on
[freq, fftResult, psdResult] = data2fftpsd(wodata(2,:),Fs);
plot(freq,fftResult);
grid on
xlabel('hz');
ylabel('magnitude')
legend('with','with out')

figure(2)
title('PSD')
clf
[freq, fftResult, psdResult] = data2fftpsd(wdata(2,:),Fs);
plot(freq,psdResult);
hold on
[freq, fftResult, psdResult] = data2fftpsd(wodata(2,:),Fs);
plot(freq,psdResult);
grid on
xlabel('hz');
ylabel('psd')
legend('with','with out')

mode = 3;

%% data selection
t = wdata(1,:);
y = wdata(2,:);

%%

if mode == 3

% test signal
% x = 0.01:0.01:10;
% y = sin(2*pi*0.3*x) + 2 * sin(2*pi*5*x);
% yinte = -cos(2*pi*0.3*x)/(2*pi*0.3) - 2 * cos(2*pi*5*x) / (2*pi*5);
% yinte2 = -sin(2*pi*0.3*x)/(2*pi*0.3)^2 -2*sin(2*pi*5*x)/(2*pi*5)^2;

% regen engine
L = length(y);
A = fftshift(fft(y));
df = Fs/L;

if ~mod(L,2)
    f = df*(-L/2:L/2-1); % n is even
else
    f = df*(-(L-1)/2:(L-1)/2); % n is odd
end

wk = 2*pi*f;
wk((wk==0)) = 1e10;

AA = zeros(1,L);
V = zeros(1,L);
D = zeros(1,L);
for i = 1:length(A)
    AA(i) = A(i) *1i*wk(i);
    V(i) = A(i) *-1i/wk(i);
    D(i) = V(i) *-1i/wk(i);
end

aa = ifft(ifftshift(AA),'symmetric');
v = ifft(ifftshift(V),'symmetric');
d = ifft(ifftshift(D),'symmetric');

v = v - v(1);
d = d - d(1);

figure(3)
clf
hold on
plot(t,y)
plot(t,v)
plot(t,d)
grid on

figure(4)
clf
hold on
plot(t,v)
plot(t,cumsum(y/Fs))
% plot(t,v-cumsum(y/Fs))

figure(5)
clf
hold on
plot(t,d)
plot(t,cumsum(cumsum(y/Fs)/Fs))
% plot(t,d-cumsum(cumsum(y/Fs)/Fs))

end