function out = findNearestIdx(data, arg)

dist = (data-arg).^2;

[~,out] = min(dist);

end