function out = getActiveVeh(tick, stuckVeh)

tick(stuckVeh==true) = inf;
[~,out] = min(tick);

end