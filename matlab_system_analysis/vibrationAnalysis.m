%% Load vibration data
temp = load('wDamper_01.mat');
wdata = temp.Low_g_Acceleration;
wdata = wdata(:,9000:60000);
temp = load('woDamper_01.mat');
wodata = temp.Low_g_Acceleration;
% wodata = wodata(:,10000:50000);

figure(10)
clf
plot(wdata(2,:))

for i = 1:3
    wdata(i+1,:) = detrend(wdata(i+1,:));
    wodata(i+1,:) = detrend(wodata(i+1,:));
end

Fs = length(wdata)/wdata(1,end);

figure(1)
title('FFT')
clf
[freq, fftResult, psdResult] = data2fftpsd(wdata(2,:),Fs);
plot(freq,fftResult);x`
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

mode = 1;


%%

x = wdata(1,:);
y = wdata(4,:);

% x = 0:0.001:10;
% y = cos(2*pi*x);
% y = x - floor(x);
% Fs = 1000;
% y_comp = sin(2*pi*x)/2/pi;
% y_comp = y.^2/2 + floor(x)/2;

fftY = fft(y);

[freq_fft, fftResult]=data2fftpsd(y,Fs);

L = length(fftResult);

% figure(3)
% plot(freq,fftResult)
% title('fft')

%%
if mode == 1

Lraw = length(x);
testSig = zeros(1,Lraw);
Sig_inte = zeros(1,Lraw);
Sig_diff = zeros(1,Lraw);
Sig_inte2 = zeros(1,Lraw);
freq = Fs*(0:length(x)-1) / Lraw;
radius = abs(fftY)/Lraw;
phase = angle(fftY);

prevDone = 0;
curDone = 0;

% % regen engine
for i = 1:Lraw
    testSig = testSig + radius(i) * sin(2*pi*freq(i)*x+phase(i)+pi/2);
    
%     curDone = floor(i/Lraw*100);
%     if (prevDone ~= curDone)
%         disp([num2str(curDone), '%']);
%         prevDone = curDone;
%     end
    curDone = floor(i/Lraw*100);
    if (prevDone ~= curDone)
        text = '[';
        for n = 1:curDone
            text = text + "#";
        end
        for n = 1:100-curDone
            text = text + "_";
        end
        text = text + "] " + num2str(curDone) + "%";
        disp(text)
        prevDone = curDone;
    end
end

% integration engine
for i = 1:Lraw
    if freq(i)~=0
        Sig_inte = Sig_inte - radius(i) / (2*pi*freq(i)) * cos(2*pi*freq(i)*x + phase(i) + pi/2);
    else
        Sig_inte = Sig_inte + radius(i) * sin(phase(i)+pi/2) * x;
    end

    curDone = floor(i/Lraw*100);
    if (prevDone ~= curDone)
        text = '[';
        for n = 1:curDone
            text = text + "#";
        end
        for n = 1:100-curDone
            text = text + "_";
        end
        text = text + "] " + num2str(curDone) + "%";
        disp(text)
        prevDone = curDone;
    end
end
Sig_inte = Sig_inte - Sig_inte(1);

% differential engine
for i = 1:Lraw
%         testSig = testSig + radius(i) * sin(2*pi*freq(i)*x+phase(i)+pi/2);
    Sig_diff = Sig_diff + radius(i) * 2*pi*freq(i) * cos(2*pi*freq(i)*x + phase(i)+pi/2);

    curDone = floor(i/Lraw*100);
    if (prevDone ~= curDone)
        text = '[';
        for n = 1:curDone
            text = text + "#";
        end
        for n = 1:100-curDone
            text = text + "_";
        end
        text = text + "] " + num2str(curDone) + "%";
        disp(text)
        prevDone = curDone;
    end
end

% double integration engine
for i = 1:Lraw
    if freq(i)~=0
        Sig_inte2 = Sig_inte2 - radius(i) / (2*pi*freq(i))^2 * sin(2*pi*freq(i)*x + phase(i) + pi/2);
    else
        Sig_inte2 = Sig_inte2 + radius(i) * sin(phase(i)+pi/2) * x;
    end    

    curDone = floor(i/Lraw*100);
    if (prevDone ~= curDone)
        text = '[';
        for n = 1:curDone
            text = text + "#";
        end
        for n = 1:100-curDone
            text = text + "_";
        end
        text = text + "] " + num2str(curDone) + "%";
        disp(text)
        prevDone = curDone;
    end
end
Sig_inte2 = Sig_inte2 - Sig_inte2(1);


end

%%

if mode == 2
L = length(fftResult);
Lraw = length(x);
prevDone = 0;
curDone = 0;
testSig = zeros(1,Lraw);
Sig_inte = zeros(1,Lraw);
Sig_diff = zeros(1,Lraw);

% regen engine
for i = 1:L

    testSig = testSig + fftResult(i) * sin(2*pi*freq_fft(i)*x);

    curDone = floor(i/L*100);
    if (prevDone ~= curDone)
        text = '[';
        for n = 1:curDone
            text = text + "#";
        end
        for n = 1:100-curDone
            text = text + "_";
        end
        text = text + "] " + num2str(curDone) + "%";
        disp(text)
        prevDone = curDone;
    end
end

% % integration engine
for i = 1:L
    if freq(i)~=0
        Sig_inte = Sig_inte - fftResult(i) / (2*pi*freq_fft(i)) * cos(2*pi*freq_fft(i)*x);
    else
    end

    curDone = floor(i/L*100);
    if (prevDone ~= curDone)
        text = '[';
        for n = 1:curDone
            text = text + "#";
        end
        for n = 1:100-curDone
            text = text + "_";
        end
        text = text + "] " + num2str(curDone) + "%";
        disp(text)
        prevDone = curDone;
    end
end
% 
% % differential engine
for i = 1:L
%     Sig_diff = Sig_diff + radius(i) * 2*pi*freq(i) * cos(2*pi*freq(i)*x + phase(i)+pi/2);
    Sig_diff = Sig_diff + fftResult(i) * 2*pi*freq_fft(i) * cos(2*pi*freq_fft(i)*x);

    curDone = floor(i/L*100);
    if (prevDone ~= curDone)
        text = '[';
        for n = 1:curDone
            text = text + "#";
        end
        for n = 1:100-curDone
            text = text + "_";
        end
        text = text + "] " + num2str(curDone) + "%";
        disp(text)
        prevDone = curDone;
    end
end
end
%%
figure(4)
clf
plot(x,Sig_inte)
hold on
% plot(x,y_comp)
plot(x,cumsum(y)/Fs)
title('integrated signal')

figure(5)
clf
plot(x,y)
title('original signal')

figure(6)
clf
plot(x,testSig)
title('regenerated signal')

figure(7)
clf
plot(x,Sig_diff)
hold on
plot(x(1:end-1),diff(y)*Fs)
title('differentiated signal')

figure(8)
clf
plot(x,Sig_inte2)
hold on
plot(x,cumsum(cumsum(y)/Fs)/Fs)
title('double integration signal')
