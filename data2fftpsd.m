function out = data2fftpsd(data)

X = data;
X = X-mean(X);
L = length(X);

Fs = 50; % sampling rate

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:L/2)/L;
plot(f(1:end),P1(1:end))
hold on

Pxx = 1/(L*Fs)*abs(Y(1:length(X)/2+1)).^2;
freq = 0:Fs/L:Fs/2;
semilogy(freq(2:end),Pxx(2:end));

grid on
legend('FFT','PSD')

end