%% Load File
gdLog = '/home/jaehan/Desktop/test_flight/gdLogCsv/gdLog_200409_102000.csv';
[data,time] = loader(gdLog);

%% Parameter Setting
range = [19530 20610];          % Specify data range
freq_range = [0 1/4];           % Specify frequency range for signal filtering [hz]

%% Main
T = time(range(1):range(2));
X = zeromean(data.rpy_0(range(1):range(2)));
L = length(X);
Fs = 50;                        % System Hz

% Original Signal
figure(1)
clf
plot(T,X)
grid on
title('Original signal','fontsize',14)
ylabel('m','fontsize',14)

% Fourier Transformed plot
figure(2)
clf
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:L/2)/L;
plot(f(1:end),P1(1:end))
title('FFT acceleration','fontsize',14)
xlabel('Hz','fontsize',14)
ylabel('m/s^2','fontsize',14)
grid on
hold on

% Power Spectral Density plot
figure(3)
clf
Pxx = 1/(L*Fs)*abs(Y(1:length(X)/2+1)).^2;
freq = 0:Fs/L:Fs/2;
semilogy(freq(2:end),Pxx(2:end));
xlabel('Hz'); ylabel('(m/s^2)^2/Hz');
title('PSD','fontsize',14)
grid on
hold on
plot([15 15],[0 0.1],'k--')

%% Signal filtering

boola = f>=freq_range(1);
boolb = f<=freq_range(2);
boolf = boola.*boolb;

P1_rev = P1.*boolf';

boolf_2 = [boolf,fliplr(boolf(1:end-1))];
Y_rev = Y.*boolf_2';
X_rev = ifft(Y_rev);

figure(1)
hold on
plot(T,real(X_rev),'r');
title('Frequency Discrimination','fontsize',14)
ylabel('m','fontsize',14)
grid on
legend('original signal','filtered')



%% animation (for fun!)
% 
% t_sample = seconds(time);
% t_sample = t_sample(1:600);
% d_sample = real(X_rev);
% d_sample = d_sample(1:600);
% 
% figure()
% plot([t_sample(1) t_sample(end)],[0 0],'k')
% hold on
% curve = animatedline;
% 
% xlim([t_sample(1) t_sample(end)]);
% ylim([-1.2 1.2])
% grid on
% 
% for i = 1:length(d_sample)
% addpoints(curve,t_sample(i),d_sample(i));
% drawnow
% pause(0.009)
% end