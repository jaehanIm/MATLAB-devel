function out = slopeQuadResponseSurf(coeff,x,y)
a = coeff(1);
b = coeff(2);
c = coeff(3);
d = coeff(4);
e = coeff(5);
f = coeff(6);

dx = 2*a*x + c*y + d;
dy = 2*b*y + c*x + e;
dz = -1;

out = [dx,dy,dz];
end