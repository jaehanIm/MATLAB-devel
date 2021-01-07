% close all

% gdLog = '/home/jaehan/Desktop/test flight/gdLogCsv/gdLog_200414_113807.csv'; % 용대뤼
% gdLog = '/home/jaehan/Desktop/test flight/gyeongju/gdLog_200121_094823.csv'; % 경주 data
% gdLog = '/home/jaehan/Desktop/test flight/gdLog_200520_153847.csv';
% gdLog = '/home/jaehan/Desktop/test flight/200715_154607/gdLog_200715_154607.csv';
% [time, out, acc, gimbal] = logReader(gdLog);

% temporary files
gdLog = '/home/jaehan/Desktop/test flight/Vehicle_Analysis/KH/201013_105028/gdLog_201013_105028.csv';

[data,time] = loader(gdLog);

% sample = [4600 8440]; % 1st lidar
% sample = [8440 9678]; % 1st hold
% sample = [9678 13170]; % 2nd lidar
% sample = [14190 19340]; % 3rd zlidar
sample = [19530 20610]; % temporary

T = time(sample(1):sample(2));
% X = zeromean(acc(sample(1):sample(2),2));
X = zeromean(data.gimbalRPY_2(sample(1):sample(2)));
% X = zeromean([0;diff(data.rpy_0(sample(1):sample(2)))]);
% X = zeromean(rot_data(1,sample(1):sample(2)))';
L = length(X);
Fs = 50;

figure(1)
clf
plot(T,X)
grid on
title('Original signal','fontsize',14)
ylabel('m','fontsize',14)


figure(2)
clf
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:L/2)/L;
plot(f(1:end),P1(1:end))
% semilogy(f,P1)
title('FFT acceleration','fontsize',14)
xlabel('Hz','fontsize',14)
ylabel('m/s^2','fontsize',14)
grid on
hold on
plot([15 15],[0 0.1],'k--')

figure(3)
clf
Pxx = 1/(L*Fs)*abs(Y(1:length(X)/2+1)).^2;
freq = 0:Fs/L:Fs/2;
% plot(freq(3:end),Pxx(3:end));
semilogy(freq(2:end),Pxx(2:end));
% plot(freq(2:end),Pxx(2:end));
xlabel('Hz'); ylabel('(m/s^2)^2/Hz');
title('PSD','fontsize',14)
grid on
hold on
plot([15 15],[0 0.1],'k--')

%% Going backacc

freq_range = [0 1/4];
boola = f>=freq_range(1);
boolb = f<=freq_range(2);
boolf = boola.*boolb;

P1_rev = P1.*boolf';

boolf_2 = [boolf,fliplr(boolf(1:end-1))];
Y_rev = Y.*boolf_2';
% Y_rev(1) = Y(1); Y_rev(end) = Y(end);
X_rev = ifft(Y_rev);

figure(1)
hold on
plot(T,abs(X_rev).*real(X_rev)./abs(real(X_rev)),'r');
title('Frequency Discrimination','fontsize',14)
ylabel('m','fontsize',14)
grid on
% xlim([T(1) T(600)])
% ylim([-1.2 1.2])
legend('original signal','filtered')



%% animation
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