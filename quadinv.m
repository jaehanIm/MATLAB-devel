function out = quadinv(f,x)

a = f(1);
b = f(2);
c = f(3);

out(1) = (-b+sqrt(b^2-4*a*(c-x)))/(2*a);
out(2) = (-b-sqrt(b^2-4*a*(c-x)))/(2*a);

end