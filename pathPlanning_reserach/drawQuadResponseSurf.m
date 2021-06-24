function out = drawQuadResponseSurf(coeff,X,Y)
a = coeff(1);
b = coeff(2);
c = coeff(3);
d = coeff(4);
e = coeff(5);
f = coeff(6);

out = a*X.^2 + b*Y.^2 + c*X.*Y + d*X + e*Y + f;
end