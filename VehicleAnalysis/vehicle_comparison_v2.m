tic

%% Data selection

addpath('./1013FT')

clear tfDataSet
tfDataSet = [];

temp = load('KH_iodelay.mat'); % 대조군
% SB = load('KH_0918'); SB = SB.tfResult; % 실험군
% SB = load('SB'); SB = SB.tfResult; % 실험군
% SB = load('KH_discrete'); SB = SB.tfResult; % 실험군
% SB = load('KH_sweep'); SB = SB.tfResult; % 실험군
tfDataSet = horzcat(tfDataSet,temp.tfResult');

temp = load('SB_default_iodelay');
tfDataSet = horzcat(tfDataSet,temp.tfResult');
% temp = load('SB_Prop_iodelay');
% tfDataSet = horzcat(tfDataSet,temp.tfResult');
% temp = load('SB_advanced_iodelay');
% tfDataSet = horzcat(tfDataSet,temp.tfResult');
% temp = load('SB_cpp_iodelay');
% tfDataSet = horzcat(tfDataSet,temp.tfResult');

Flag = ["Roll","v","Y","Pitch","u","X","Yaw","w","Alt"];
dataFlag = ["KH","SB default","SB GainSet","SB advanced","SB robustness"];
% dataFlag = ["KH","SB GainSet"];
% dataFlag = ["SB default","SB GainSet","SB advanced","SB robustness"];

%% Analysis

n = size(tfDataSet,2);

BW = zeros(n,9);
stepResult{n,9} = [];
stepSse = zeros(n,9);
stepDamp = zeros(n,9);

% figure(4)
% clf
% DRB = zeros(9,2,n);
% DRP = zeros(9,2,n);
% for i = 1:9
%     subplot(3,3,i)
%     hold on
%     for j = 1:n
%         curTf = tf(tfDataSet{i,j}.Num,tfDataSet{i,j}.Den);
%         G = -curTf/(curTf-1);
%         S = 1/(1+G);
%         bode(S);
%         DRB(i,:,j) = bandwidth(S);
%         DRP(i,:,j) = mag2db(getPeakGain(S));
%         drawnow;
%     end
%     grid on
%     title("DRB "+Flag(i))
% end
% legend(dataFlag)
% disp('c0')

figure(5)
clf
for i = 1:9
    subplot(3,3,i)
    hold on
    for j = 1:n
        curTf = tf(tfDataSet{i,j}.Num,tfDataSet{i,j}.Den,'iodelay',tfDataSet{i,j}.Delay);
        y = step(curTf);
        step(curTf);
        title("Step "+Flag(i))
        BW(j,i) = bandwidth(curTf);
        stepResult{j,i} = stepinfo(curTf);
        stepSse(j,i) = 1 - y(end);
        [~,temp] = damp(curTf);
        stepDamp(j,i) = temp(1);
        drawnow
    end
    grid on
end
legend(dataFlag)
disp('c1')

figure(6)
clf
GMPM = zeros(9,2,n);
for i = 1:9
    subplot(3,3,i)
    hold on
    for j = 1:n
        [Gm,Pm] = margin(tf(tfDataSet{i,j}.Num,tfDataSet{i,j}.Den));
        margin(tf(tfDataSet{i,j}.Num,tfDataSet{i,j}.Den));
        GMPM(i,:,j) = [mag2db(Gm),Pm];
        drawnow;
    end
    grid on
    xlim([1 31.5])
    title("Bode "+Flag(i))
end
legend(dataFlag)
disp('c2')
%%mix code
% figure(7)
% clf
% for i = 1:3
%     subplot(3,1,i)
%     hold on
%     grid on
%     step(tf(KH_mix{i}.Num,KH_mix{i}.Den))
%     step(tf(SB_mix{i}.Num,SB_mix{i}.Den))
%     title("Step "+FlagMix(i))
%     BW_KH_mix(i) = bandwidth(tf(KH_mix{i}.Num,KH_mix{i}.Den));
%     BW_SB_mix(i) = bandwidth(tf(SB_mix{i}.Num,SB_mix{i}.Den));
%     stepResult{1,i+9} = stepinfo(tf(KH_mix{i}.Num,KH_mix{i}.Den));
%     stepResult{2,i+9} = stepinfo(tf(SB_mix{i}.Num,SB_mix{i}.Den));
%     
%     grid on
% end
% disp('c3')
% figure(8)
% clf
% GMPM_KH_mix = zeros(3,4);
% GMPM_SB_mix = zeros(3,4);
% for i = 1:3
%     subplot(3,1,i)
%     hold on
%     grid on
%     [Gm,Pm,Wcg,Wcp]=margin(tf(KH_mix{i}.Num,KH_mix{i}.Den));
%     GMPM_KH_mix(i,:) = [mag2db(Gm),Pm,Wcg,Wcp];
%     margin(tf(KH_mix{i}.Num,KH_mix{i}.Den));
%     [Gm,Pm,Wcg,Wcp]=margin(tf(SB_mix{i}.Num,SB_mix{i}.Den));
%     GMPM_SB_mix(i,:) = [mag2db(Gm),Pm,Wcg,Wcp];
%     margin(tf(SB_mix{i}.Num,SB_mix{i}.Den));
%     title("Bode "+FlagMix(i))
%     grid on
% end
% disp('c4')

stepRise = zeros(n,9);
stepSettle = zeros(n,9);
stepOver = zeros(n,9);
for i = 1:9
    for j = 1:n
        stepRise(j,i) = stepResult{j,i}.RiseTime;
        stepSettle(j,i) = stepResult{j,i}.SettlingTime;
        stepOver(j,i) = stepResult{j,i}.Overshoot;
    end
end


figure(9)
clf
subplot(6,1,1)
% bar([[BW_KH,BW_KH_mix]'/pi/2,[BW_SB,BW_SB_mix]'/pi/2])
bar(BW')
xticklabels(Flag)
grid on
title('BandWidth')
ylabel('hz')
legend(dataFlag)

subplot(6,1,2)
% bar([[GMPM_KH(:,1);GMPM_KH_mix(:,1)],[GMPM_SB(:,1);GMPM_SB_mix(:,1)]])
bar(reshape(GMPM(:,1,:),9,n))
xticklabels(Flag)
grid on
title('Gain Margin')
ylabel('db')

subplot(6,1,3)
bar(reshape(GMPM(:,2,:),9,n))
xticklabels(Flag)
grid on
title('Phase Margin')
ylabel('deg')

subplot(6,1,4)
bar(stepRise')
xticklabels(Flag)
grid on
title('Rise Time')
ylabel('s')

subplot(6,1,5)
bar(stepSettle')
xticklabels(Flag)
grid on
title('Settling Time')
ylabel('s')

subplot(6,1,6)
bar(abs(stepSse'))
xticklabels(Flag)
grid on
title('Overshoot')
ylabel('%')
toc


figure(10)
clf
subplot(6,1,1)
temp = horzcat(BW')/2/pi;
mtemp = max(temp,[],2);
temp = temp./mtemp;
bar(temp)
xticklabels(Flag)
grid on
title('BandWidth')
ylabel('hz')
legend(dataFlag)

subplot(6,1,2)
temp = reshape(GMPM(:,1,:),9,n);
mtemp = max(temp,[],2);
temp = temp./mtemp;
bar(temp)
xticklabels(Flag)
grid on
title('Gain Margin')
ylabel('db')

subplot(6,1,3)
temp = reshape(GMPM(:,2,:),9,n);
mtemp = max(temp,[],2);
temp = temp./mtemp;
bar(temp)
xticklabels(Flag)
grid on
title('Phase Margin')
ylabel('deg')

subplot(6,1,4)
temp = stepRise';
mtemp = max(temp,[],2);
temp = temp./mtemp;
bar(temp)
xticklabels(Flag)
grid on
title('Rise Time')
ylabel('s')

subplot(6,1,5)
temp = stepSettle';
mtemp = max(temp,[],2);
temp = temp./mtemp;
bar(temp)
xticklabels(Flag)
grid on
title('Settling Time')
ylabel('s')

subplot(6,1,6)
temp = abs(stepDamp');
mtemp = 1;
temp = temp./mtemp;
bar(temp)
xticklabels(Flag)
grid on
title('Damping ratio')
ylabel('%')

figure(11)
clf
for i = 1:9
    subplot(3,3,i)
    hold on
    for j = 1:n
        curTf = tf(tfDataSet{i,j}.Num,tfDataSet{i,j}.Den);
        rlocus(curTf)
        drawnow;
    end
    grid on
    title("Root Locus "+Flag(i))
end
legend(dataFlag)
disp('c4')
