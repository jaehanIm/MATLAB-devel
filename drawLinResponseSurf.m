function out = drawLinResponseSurf(coeff,X,Y)
a = coeff(1);
b = coeff(2);
c = coeff(3);

out = a*X + b*Y + c;
end