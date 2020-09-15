KH = load('KH'); KH = KH.tfResult;
SB = load('SB'); SB = SB.tfResult;
KH_mix = load('KH_mix'); KH_mix = KH_mix.tfResult;
SB_mix = load('SB_mix'); SB_mix = SB_mix.tfResult;
figure(5)
clf
for i = 1:9
    subplot(3,3,i)
    hold on
    grid on
    step(tf(KH{i}.Num,KH{i}.Den))
    step(tf(SB{i}.Num,SB{i}.Den))
end
disp('c1')

figure(6)
clf
for i = 1:9
    subplot(3,3,i)
    hold on
    grid on
    margin(tf(KH{i}.Num,KH{i}.Den))
    margin(tf(SB{i}.Num,SB{i}.Den))
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
end
disp('c3')
figure(8)
clf
for i = 1:3
    subplot(3,1,i)
    hold on
    grid on
    margin(tf(KH_mix{i}.Num,KH_mix{i}.Den))
    margin(tf(SB_mix{i}.Num,SB_mix{i}.Den))
end
disp('c4')