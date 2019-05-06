function dz  = dynFunc(~,z,mu)
dz(2,1) = 0;
dz(1) = z(2);
dz(2) = -mu*(z(1)^2-1)*z(2)-z(1);
end