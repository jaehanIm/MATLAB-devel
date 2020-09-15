% close all
d2r = pi/180; r2d = 1/d2r; g = 9.81; % d2r & 중력상수

%게인값 설정
kp_pos = 0.85; kd_pos = 0.675;
kp_vel = 12 * d2r; kd_vel = 0;
w = 0.01:0.01:6.28*100; % test 대역 (0~10hz)

% pitch angle response parameter (system ID.에서 나온 결과)
wn = sqrt(94.77);
k = 82.85/wn^2;
zeta = 9.124/2/wn; 
% wn = sqrt(636.6);
% k = 626/wn^2;
% zeta = 321.5/2/wn; 

gain = wn^2./sqrt((wn^2-w.^2).^2+4.*zeta^2.*wn^2.*w.^2)*k;  % pitch angle response gain
f = w/2/pi;                                                 % rad/s to hz
gain_pos = sqrt(kp_pos^2+(w*kd_pos).^2);                    % position PID controller의 gain
gain_vel = sqrt(kp_vel^2+(w*kd_vel).^2);                    % velocity PID controller의 gain
gain_int = 1./w;                                            % 적분기 gain
gain_d = w;
% Total controller(postion response)의 gain
gain_pa = (gain_pos.*gain_vel.*gain.*(gain_int.^2)*g)./(1+gain_vel.*gain.*gain_int*g+gain_pos.*gain_vel.*gain.*(gain_int.^2)*g).*gain_d.^2/g;
gain_pp = (gain_pos.*gain_vel.*gain.*(gain_int.^2)*g)./(1+gain_vel.*gain.*gain_int*g+gain_pos.*gain_vel.*gain.*(gain_int.^2)*g);

% A->A
figure(1)
clf
semilogx(f,20*log(gain))
grid on
title('Gain plot (A_c_m_d->A)')
xlabel('hz')
ylabel('db')

% P->A
figure(2)
clf
semilogx(f,20*log(gain_pa))
grid on
title('Gain plot (P_c_m_d->A)')
xlabel('hz')
ylabel('db')
f(1091)

% P->P
figure(3)
clf
semilogx(f,20*log(gain_pp))
grid on
title('Gain plot (P_c_m_d->P)')
xlabel('hz')
ylabel('db')

% A->P
figure(4)
clf
semilogx(f,20*log(1./gain_pa.*gain))
grid on
title('Gain plot (A_c_m_d->P)')
xlabel('hz')
ylabel('db')

% Anything
figure(5)
clf
semilogx(f,20*log(1./gain_pa.*gain_pp)) % A->P
% semilogx(f,20*log(1./gain)) % A->Acmd
grid on
title('Gain plot (A->P)')
xlabel('hz')
ylabel('db')

%% Peak frequency response along with cutoff frequency

%%% This part is legacy %%%
% wn_p = sqrt(1.066);
% k_p = 1.105/wn_p^2;
% zeta_p = 2.265/2/wn_p; % position command
%%% ------------------- %%%

% 변수 변경 (pitch angle response)
wn_a = wn;
k_a = k;
zeta_a = zeta; 

f_pa = 1.7364 * 2 * pi;
f_pp = 0.2 * 2 * pi;
% Angle response at 1.682hz (피치각 반응의 peak freq.임)
gain_a_pa = wn_a^2./sqrt((wn_a^2-f_pa^2)^2+4*zeta_a^2*wn_a^2*f_pa^2)*k_a;
% Angle response at 1.151hz (위치 반응의 peak freq.임)
gain_a_pp = wn_a^2./sqrt((wn_a^2-f_pp^2)^2+4*zeta_a^2*wn_a^2*f_pp^2)*k_a;

% Position response at 0.02 hz (위치 반응의 peak freq.임)
gain_pos = sqrt(kp_pos^2+(f_pa*kd_pos).^2);
gain_vel = sqrt(kp_vel^2+(f_pa*kd_vel).^2);
gain_int = 1./f_pa;
gain_d = f_pa;
gain_pa = (gain_pos.*gain_vel.*gain_a_pa.*(gain_int.^2)*g)./(1+gain_vel.*gain_a_pa.*gain_int*g+gain_pos.*gain_vel.*gain_a_pa.*(gain_int.^2)*g).*gain_d.^2/g;
gain_pos = sqrt(kp_pos^2+(f_pp*kd_pos).^2);
gain_vel = sqrt(kp_vel^2+(f_pp*kd_vel).^2);
gain_int = 1./f_pp;
gain_d = f_pp;
gain_pp = (gain_pos.*gain_vel.*gain_a_pp.*(gain_int.^2)*g)./(1+gain_vel.*gain_a_pp.*gain_int*g+gain_pos.*gain_vel.*gain_a_pp.*(gain_int).^2*g);
% gain_fb = 1./gain_fb;
% cutoff hz test range (0~4hz)
cutoff_hz = 0.1:0.01:4;
result = zeros(length(cutoff_hz),3);
n = 1;
for cutoff_hz = 0.1:0.01:5
    wc = cutoff_hz*2*pi;
    filter_low = 1./sqrt(1+(f_pp/wc).^2); % 1st order filter gain @ f_pp
    filter_high = 1./sqrt(1+(f_pa/wc).^2); % 1st order filter gain @ f_pa
    gain_low = filter_low * gain_pp;
    gain_high = filter_high * gain_pa;
    result(n,:) = [cutoff_hz,gain_low,gain_high];
    n = n+1;
end
% 
figure(6)
% clf
hold on
plot(result(:,1),normalize(result(:,2),'range'),'b')
plot(result(:,1),normalize(result(:,3),'range'),'r')
plot([0.3 0.3],[0 1.2],'k')
plot([0.5 0.5],[0 1.2],'k')
plot([0.7 0.7],[0 1.2],'k')
xlabel('Cutoff freq.')
ylabel('Relative gain')
legend('P->P response','P->A response')
title('Normalized gain')
% % % 

% Cutoff freq to gain @ peak frequencies plot
figure(7)
% clf
hold on
plot(result(:,1),result(:,2),'b')
plot(result(:,1),result(:,3),'r')
plot([0.3 0.3],[0 1.2],'k')
plot([0.5 0.5],[0 1.2],'k')
plot([0.7 0.7],[0 1.2],'k')
legend('P->P response @0.2hz','P->A response @0.73hz')
grid on
title('Most sensitive frequency response','fontsize',14)
xlabel('Cutoff frequency')
ylabel('Gain')
% set(gca, 'YScale', 'log')

figure(8)
% clf
hold on
plot(result(:,1),abs(normalize(result(:,2),'range')-normalize(result(:,3),'range')),'b')
% plot(result(:,1),abs(result(:,2)-result(:,3)))
plot([0.3 0.3],[0 .8],'k')
plot([0.5 0.5],[0 .8],'k')
plot([0.7 0.7],[0 .8],'k')
grid on
title('Gain difference (Normalized)')
xlabel('Cutoff freq.')
ylabel('gain')

figure(9)
% clf
hold on
% plot(result(:,1),abs(normalize(result(:,2),'range')-normalize(result(:,3),'range')),'b')
plot(result(:,1),abs(result(:,2)-result(end,2)+result(end,3)-result(:,3)),'b')
plot([0.3 0.3],[0 .8],'k')
plot([0.5 0.5],[0 .8],'k')
plot([0.7 0.7],[0 .8],'k')
grid on
title('Gain difference (Normalized)')
xlabel('Cutoff freq.')
ylabel('gain')

figure(10)
% clf
hold on
plot(result(:,1),result(:,2)-result(end,2)+result(end,3),'b')
plot(result(:,1),result(:,3),'r')
plot([0.3 0.3],[0 1.1],'k')
plot([0.5 0.5],[0 1.1],'k')
plot([0.7 0.7],[0 1.1],'k')
legend('P->P response','P->A response')
grid on
title('Most sensitive frequency response','fontsize',14)
xlabel('Cutoff frequency')
ylabel('Gain')
% set(gca, 'YScale', 'log')