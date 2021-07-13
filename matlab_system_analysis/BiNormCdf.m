function out = BiNormCdf(X,Y,mu,sigma,rangeX,rangeY)

pdf = BiNorm(X,Y,mu,sigma);

% grid size recognition
yStep = Y(2,1)-Y(1,1);
xStep = X(1,2)-X(1,1);
gridArea = xStep*yStep;

% find center

[~,center_X_index] = min(abs(X-rangeX(1)),[],2); center_X_index = center_X_index(1);
[~,center_Y_index] = min(abs(Y-rangeY(1)),[],1); center_Y_index = center_Y_index(1);

% find X, Y fix satisfying rangeX, rangeY
[~,range_X_index] = min(abs(X-rangeX(2)),[],2); range_X_index = range_X_index(1);
[~,range_Y_index] = min(abs(Y-rangeY(2)),[],1); range_Y_index = range_Y_index(1);

% integrate pdf and finish
xi = min(center_X_index,range_X_index);
xf = max(center_X_index,range_X_index);
yi = min(center_Y_index,range_Y_index);
yf = max(center_Y_index,range_Y_index);

range_pdf = pdf(yi:yf,xi:xf)*gridArea;
out = sum(range_pdf,'all');

end