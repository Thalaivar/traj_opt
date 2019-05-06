function dz  = dynFunc(~,z,a,b)

dz(2,1) = 0;
dz(1) = (-z(1) + a*z(2) + (z(1)^2)*z(2));
dz(2) = (b - a*z(2) - (z(1)^2)*z(2));

end