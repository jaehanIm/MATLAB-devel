function out = resample(in,step)
n = floor(length(in)/step);

out = zeros(n,1);
for i = 1:n
    out(i) = in(1+step*(i-1));
end

end