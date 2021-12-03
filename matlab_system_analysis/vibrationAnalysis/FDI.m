function v = FDI(y,Fs)

% regen engine
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

V = zeros(1,L);
for i = 1:length(A)
    V(i) = A(i) *-1i/wk(i);
end

v = ifft(ifftshift(V),'symmetric');
v = v - v(1);
v = v';

end