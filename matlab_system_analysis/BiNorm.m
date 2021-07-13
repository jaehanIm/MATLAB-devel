%BiNorm.m
% Plots the Bivariate Normal Density Function for Various Parameters Values
% author: David Dobor 

function [pdf] = BiNorm(x,y,mu,sigma)

const = (1/sqrt(2*pi))^2;
const = const/sqrt(det(sigma));
temp = [x(:)-mu(1) y(:)-mu(2)];
pdf = const*exp(-0.5*diag(temp/(sigma)*temp'));
pdf = reshape(pdf,size(x));

end