KH = load('KH'); KH = KH.tfResult; % 대조군
SB = load('KH_0918'); SB = SB.tfResult; % 실험군
% SB = load('SB'); SB = SB.tfResult; % 실험군
% SB = load('KH_discrete'); SB = SB.tfResult; % 실험군
% SB = load('KH_sweep'); SB = SB.tfResult; % 실험군

KH_mix = load('KH_mix'); KH_mix = KH_mix.tfResult;
SB_mix = load('KH_mix_0918'); SB_mix = SB_mix.tfResult;
Flag = ["Roll","v","Y","Pitch","u","X","Yaw","w","Alt"];
FlagMix = ["P->A","P->V","P->P"];

BW_KH = zeros(1,9); BW_SB = zeros(1,9);
BW_KH_mix = zeros(1,3); BW_SB_mix = zeros(1,3);
stepResult{2,12} = [];
stepSse = zeros(2,12);

figure(5)
clf
for i = 1:9
    subplot(3,3,i)
    hold on
    grid on
    y1 = step(tf(KH{i}.Num,KH{i}.Den));
    y2 = step(tf(SB{i}.Num,SB{i}.Den));
    step(tf(KH{i}.Num,KH{i}.Den));
    step(tf(SB{i}.Num,SB{i}.Den));
    title("Step "+Flag(i))
    BW_KH(i) = bandwidth(tf(KH{i}.Num,KH{i}.Den));
    BW_SB(i) = bandwidth(tf(SB{i}.Num,SB{i}.Den));
    stepResult{1,i} = stepinfo(tf(KH{i}.Num,KH{i}.Den));
    stepResult{2,i} = stepinfo(tf(SB{i}.Num,SB{i}.Den));
    stepSse(1,i) = 1 - y1(end);
    stepSse(2,i) = 1 - y2(end);
end
disp('c1')

figure(6)
clf
GMPM_KH = zeros(9,4);
GMPM_SB = zeros(9,4);
for i = 1:9
    subplot(3,3,i)
    hold on
    [Gm,Pm,Wcg,Wcp] = margin(tf(KH{i}.Num,KH{i}.Den));
    margin(tf(KH{i}.Num,KH{i}.Den));
    GMPM_KH(i,:) = [Gm,Pm,Wcg,Wcp];
    [Gm,Pm,Wcg,Wcp] = margin(tf(SB{i}.Num,SB{i}.Den));
    margin(tf(SB{i}.Num,SB{i}.Den));
    GMPM_SB(i,:) = [Gm,Pm,Wcg,Wcp];
    grid on
    title("Bode "+Flag(i))
end
disp('c2')
%%mix code
figure(7)
clf
for i = 1:3
    subplot(3,1,i)
    hold on
    grid on
    step(tf(KH_mix{i}.Num,KH_mix{i}.Den))
    step(tf(SB_mix{i}.Num,SB_mix{i}.Den))
    title("Step "+FlagMix(i))
    BW_KH_mix(i) = bandwidth(tf(KH_mix{i}.Num,KH_mix{i}.Den));
    BW_SB_mix(i) = bandwidth(tf(SB_mix{i}.Num,SB_mix{i}.Den));
    stepResult{1,i+9} = stepinfo(tf(KH_mix{i}.Num,KH_mix{i}.Den));
    stepResult{2,i+9} = stepinfo(tf(SB_mix{i}.Num,SB_mix{i}.Den));
end
disp('c3')
figure(8)
clf
GMPM_KH_mix = zeros(3,4);
GMPM_SB_mix = zeros(3,4);
for i = 1:3
    subplot(3,1,i)
    hold on
    grid on
    [Gm,Pm,Wcg,Wcp]=margin(tf(KH_mix{i}.Num,KH_mix{i}.Den));
    GMPM_KH_mix(i,:) = [Gm,Pm,Wcg,Wcp];
    margin(tf(KH_mix{i}.Num,KH_mix{i}.Den));
    [Gm,Pm,Wcg,Wcp]=margin(tf(SB_mix{i}.Num,SB_mix{i}.Den));
    GMPM_SB_mix(i,:) = [Gm,Pm,Wcg,Wcp];
    margin(tf(SB_mix{i}.Num,SB_mix{i}.Den));
    title("Bode "+FlagMix(i))
    grid on
end
disp('c4')

stepRise = zeros(2,12);
stepSettle = zeros(2,12);
stepOver = zeros(2,12);
for i = 1:12
    stepRise(1,i) = stepResult{1,i}.RiseTime;
    stepRise(2,i) = stepResult{2,i}.RiseTime;
    stepSettle(1,i) = stepResult{1,i}.SettlingTime;
    stepSettle(2,i) = stepResult{2,i}.SettlingTime;
    stepOver(1,i) = stepResult{1,i}.Overshoot;
    stepOVer(2,i) = stepResult{1,i}.Overshoot;
end


figure(9)
clf
subplot(6,1,1)
bar([[BW_KH,BW_KH_mix]'/pi/2,[BW_SB,BW_SB_mix]'/pi/2])
xticklabels([Flag,FlagMix])
grid on
title('BandWidth')
ylabel('hz')
legend('KH','SB')

subplot(6,1,2)
bar([[GMPM_KH(:,1);GMPM_KH_mix(:,1)],[GMPM_SB(:,1);GMPM_SB_mix(:,1)]])
xticklabels([Flag,FlagMix])
grid on
title('Gain Margin')
ylabel('db')

subplot(6,1,3)
bar([[GMPM_KH(:,2);GMPM_KH_mix(:,2)],[GMPM_SB(:,2);GMPM_SB_mix(:,2)]])
xticklabels([Flag,FlagMix])
grid on
title('Phase Margin')
ylabel('deg')

subplot(6,1,4)
bar([stepRise(1,:)',stepRise(2,:)'])
xticklabels([Flag,FlagMix])
grid on
title('Rise Time')
ylabel('deg')

subplot(6,1,5)
bar([stepSettle(1,:)',stepSettle(2,:)'])
xticklabels([Flag,FlagMix])
grid on
title('Settling Time')
ylabel('deg')

subplot(6,1,6)
bar([abs(stepSse(1,:)'),abs(stepSse(2,:)')])
xticklabels([Flag,FlagMix])
grid on
title('SSE')
ylabel('deg')