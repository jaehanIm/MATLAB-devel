w = 0.1:0.01:6.28*10;

% wn = sqrt(94.77);
% k = 82.85/wn^2;
% zeta = 9.124/2/wn; % pitch angle command

wn = sqrt(1.066);
k = 1.105/wn^2;
zeta = 2.265/2/wn; % position command

% wn = sqrt(5614);
% k = 5714/wn^2;
% zeta = 1561/2/wn; % position command

wnn = wn/2/pi
wd = sqrt(1-zeta^2)*wn/2/pi
gain = wn^2./sqrt((wn^2-w.^2).^2+4.*zeta^2.*wn^2.*w.^2)*k;
phase = -atan(2*zeta*wn.*w./(wn^2-w.^2));
f = w/2/pi;
figure(1)
semilogx(f,20*log(gain/gain(1)))
% semilogx(f,gain)
grid on
title('Gain plot (pitch control)')
xlabel('hz')
ylabel('db')
figure(2)
semilogx(f,phase*180/pi)
title('phase plot (pitch control)')
xlabel('hz')
ylabel('deg')
grid on

%% filter
w = 0.1:0.02:6.28*5;
wc = 30*2*pi;
gain = 1./sqrt(1+(w/wc).^2);
figure(3)
semilogx(w,gain)
grid on
title('filter')

%% iteration
wn_p = sqrt(1.066);
k_p = 1.105/wn_p^2;
zeta_p = 2.265/2/wn_p; % position command

wn_a = sqrt(94.77);
k_a = 82.85/wn_a^2;
zeta_a = 9.124/2/wn_a; % pitch angle command

gain_p = wn_p^2./sqrt((wn_p^2-0.064^2)^2+4*zeta_p^2*wn_p^2*0.064^2)*k_p;
gain_a = wn_a^2./sqrt((wn_a^2-1.135^2)^2+4*zeta_a^2*wn_a^2*1.135^2)*k_a;

hz = 0.1:0.01:4;
result = zeros(length(hz),3);
n = 1;
for hz = 0.1:0.01:4
    wc = hz*2*pi;
    filter_low = 1./sqrt(1+(0.064/wc).^2);
    filter_high = 1./sqrt(1+(1.178/wc).^2);
    gain_low = filter_low * gain_p;
    gain_high = filter_high * gain_a;
    result(n,:) = [hz,gain_low,gain_high];
    n = n+1;
end

figure(4)
clf
hold on
plot(result(:,1),normalize(result(:,2),'range'))
plot(result(:,1),normalize(result(:,3),'range'))
plot([0.5 0.5],[0 1.2],'k')
plot([0.7 0.7],[0 1.2],'k')
legend('position response','pitch angle response')

figure(5)
clf
hold on
plot(result(:,1),result(:,2)/0.9941)
plot(result(:,1),result(:,3)*1.143)
plot([0.5 0.5],[0.5 1.1],'k')
plot([0.7 0.7],[0.5 1.1],'k')
legend('position response','pitch angle response')