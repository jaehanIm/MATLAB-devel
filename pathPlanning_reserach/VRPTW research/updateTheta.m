function theta = updateTheta(N, C, currentNode, timeWindow, tick, lambda)

theta = zeros(1,N);
for j = 1:N
    startSlack = max(tick + C(currentNode,j),timeWindow(j,1));
    endSlack = timeWindow(j,2) - tick;
    theta(j) = 1/max(1,(lambda * startSlack * endSlack));
end

end