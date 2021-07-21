function [freq, fftResult, psdResult] = data2fftpsd(data, Fs)

X = data;
X = X-mean(X);
L = length(X);

Y = fft(X);
P2 = abs(Y/L);
fftResult = P2(1:L/2+1);

freq = Fs*(0:L/2)/L;
freq = freq';
fftResult(2:end-1) = 2*fftResult(2:end-1);
psdResult = 1/(L*Fs)*abs(Y(1:length(X)/2+1)).^2;


% figure(5)
% clf
% plot(freq(1:end),fftResult(1:end));
% hold on
% semilogy(freq(2:end),psdResult(2:end));
% xlim([0.1 inf])
% grid on
% legend('FFT','PSD')

end