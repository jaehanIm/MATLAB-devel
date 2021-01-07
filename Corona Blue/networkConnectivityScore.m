function score = networkConnectivityScore(agentHeight,networkHeight)
n = length(agentHeight);
score = sum(networkHeight-abs(networkHeight-agentHeight))/n/networkHeight;
end