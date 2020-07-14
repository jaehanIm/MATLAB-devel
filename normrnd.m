% Generates normally distributed random number
% Define [mu, sigma, output size (ex. [2 5] - optional)]

function out = normrnd(varargin)

mu = varargin{1};
sigma = varargin{2};
if nargin == 3
    mat = varargin{3};
else
    mat = [1,1];
end

out = zeros(mat(1),mat(2));

for i = 1:mat(1)
    for j = 1:mat(2)
        seed = rand(1);
        quantile = sqrt(2)*erfinv(2*seed-1);
        out(i,j) = mu + sigma * quantile;
    end
end

end