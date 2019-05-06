function dz  = dynFunc(~,z,A)
a = A(1); b = A(2);
dz(2,1) = 0;
dz(1) = -z(1) + a*z(2) + (z(1)^2)*z(2);
dz(2) = b - a*z(2) - (z(1)^2)*z(2);
end