% photo = readtable('/home/jaehan/log/211104/df_xmp.csv');
photo = readtable('/home/jaehan/Downloads/DSLR_quality_intersect.csv');
data = readtable('/home/jaehan/log/211104/gdLog_211104_123103.csv');

% Time processing
fTimeGlobal_Date = datetime(data.rosTime,'ConvertFrom','posixtime','TimeZone','Asia/Tokyo');
dateStart = fTimeGlobal_Date(1);

t = fTimeGlobal_Date - dateStart;
photoTime = photo.xmpTime;
photoTime = datetime(photoTime,'TimeZone','Asia/Tokyo');
photoTime = photoTime - dateStart;

t = seconds(t);
photoTime = seconds(photoTime);

%% Angle preprocessing
for i = 1:length(t)-1
    if abs(data.rpy_deg_2(i) - data.rpy_deg_2(i+1)) > 180
        adjustment = data.rpy_deg_2(i+1) - data.rpy_deg_2(i);
        data.rpy_deg_2(i+1:end) = data.rpy_deg_2(i+1:end) - adjustment; 
    end
end


%% Main
% y = sqrt(data.gimbalRpy_deg_1.^2 + data.rpy_deg_2.^2);
y = data.gimbalRpy_deg_1;
Fs = 50;

L = length(y);
A = fftshift(fft(y));
df = Fs/L;

if ~mod(L,2)
    f = df*(-L/2:L/2-1); % n is even
else
    f = df*(-(L-1)/2:(L-1)/2); % n is odd
end

wk = 2*pi*f;
wk((wk==0)) = 1e10;

% frequency filtering
% high f filter
A_orig = A;
% A(1:25000) = 0;
% A(end-25000:end) = 0;

% lowf filter
% A(L/2-3333:L/2+3333) = 0;

V = zeros(1,L);
D = zeros(1,L);
AA = zeros(1,L);
for i = 1:length(A)
    V(i) = A(i) *-1i/wk(i);
    D(i) = V(i) *-1i/wk(i);
    AA(i) = A(i) *1i*wk(i);
end

aa = ifft(ifftshift(AA),'symmetric');
v = ifft(ifftshift(V),'symmetric');
d = ifft(ifftshift(D),'symmetric');

v = v - v(1);
d = d - d(1);
% aa = aa - aa(1);

%% plot

figure(1)
clf
hold on
grid on
plot(t,aa)
plot(t,y)
plot([t(1) t(end)], [-23 -23],'k--','LineWidth',3)
plot([t(1) t(end)], [23 23],'k--','LineWidth',3)
stem(photoTime,photo.imgQuality*10,'LineWidth',5)
plot(t,data.fcMcMode*30,'m','LineWidth',1)
% xlim([110 155])
% ylim([-80 80])

figure(2)
clf 
hold on
grid on
plot(t,y)
plot(t,aa)
plot([t(1) t(end)], [-23 -23],'k--','LineWidth',3)
plot([t(1) t(end)], [23 23],'k--','LineWidth',3)
