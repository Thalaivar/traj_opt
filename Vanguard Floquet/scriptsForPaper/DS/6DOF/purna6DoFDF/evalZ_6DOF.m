function Z = evalZ_6DOF(Y,dY,wind,obj)

Wx = wind(1);
Wxz = wind(2);

phi = Y(1); thet = Y(2); psi = Y(3);
% x = Y(4); y = Y(5); z = Y(6);

dphi = dY(1,1); dthet = dY(1,2); dpsi = dY(1,3); dx = dY(1,4); dy = dY(1,5); dz = dY(1,6);
ddphi = dY(2,1); ddthet = dY(2,2); ddpsi = dY(2,3); ddx = dY(2,4); ddy = dY(2,5); ddz = dY(2,6);

Cphi = cos(phi); Sphi = sin(phi);
Cthet = cos(thet); Sthet = sin(thet);
Cpsi = cos(psi); Spsi = sin(psi);

p = dphi - dpsi*Sthet;
q = dthet*Cphi + dpsi*Cthet*Sphi;
r = -dthet*Sphi + dpsi*Cthet*Cphi;

dpdt = ddphi - dthet*dpsi*Cthet - ddpsi*Sthet;
dqdt = ddthet*Cphi - dthet*dphi*Sphi + ddpsi*Cthet*Sphi - dpsi*dthet*Sthet*Sphi + dpsi*dphi*Cthet*Cphi;
drdt = -ddthet*Sphi - dthet*dphi*Cphi + ddpsi*Cthet*Cphi - dpsi*dthet*Sthet*Cphi - dpsi*dphi*Cthet*Sphi;

u = dx*Cpsi*Cthet + dy*Spsi*Cthet - dz*Sthet;
v = dx*(Cpsi*Sthet*Sphi-Spsi*Cphi) + dy*(Spsi*Sthet*Sphi+Cpsi*Cphi) + dz*Cthet*Sphi;
w = dx*(Cpsi*Sthet*Cphi+Spsi*Sphi) + dy*(Spsi*Sthet*Cphi-Cpsi*Sphi) + dz*Cthet*Cphi;

u = u - Wx*Cpsi*Cthet;
v = v - Wx*(Cpsi*Sthet*Sphi - Spsi*Cphi);
w = w - Wx*(Cpsi*Sthet*Cphi + Spsi*Sphi);

du = ddx*Cpsi*Cthet - dx*dpsi*Spsi*Cthet - dx*dthet*Cpsi*Sthet + ddy*Spsi*Cthet + dy*dpsi*Cpsi*Cthet - dy*dthet*Spsi*Sthet - ddz*Sthet - dz*dthet*Cthet;
dv = ddx*(Cpsi*Sthet*Sphi - Spsi*Cphi) + dx*(-dpsi*Spsi*Sthet*Sphi + dthet*Cpsi*Cthet*Sphi + dphi*Cpsi*Sthet*Cphi - dpsi*Cpsi*Cphi + dphi*Spsi*Sphi)+...
    ddy*(Spsi*Sthet*Sphi + Cpsi*Cphi) + dy*(dpsi*Cpsi*Sthet*Sphi + dthet*Spsi*Cthet*Sphi + dphi*Spsi*Sthet*Cphi - dpsi*Spsi*Cphi - dphi*Cpsi*Sphi)+...
    ddz*Cthet*Sphi - dz*dthet*Sthet*Sphi + dz*dphi*Cthet*Cphi;
dw = ddx*(Cpsi*Sthet*Cphi + Spsi*Sphi) + dx*(-dpsi*Spsi*Sthet*Cphi + dthet*Cpsi*Cthet*Cphi - dphi*Cpsi*Sthet*Sphi + dpsi*Cpsi*Sphi + dphi*Spsi*Cphi)+...
    ddy*(Spsi*Sthet*Cphi - Cpsi*Sphi) + dy*(dpsi*Cpsi*Sthet*Cphi + dthet*Spsi*Cthet*Cphi - dphi*Spsi*Sthet*Sphi + dpsi*Spsi*Sphi - dphi*Cpsi*Cphi)+...
    ddz*Cthet*Cphi - dz*dthet*Sthet*Cphi - dz*dphi*Cthet*Sphi;

du = du - ( -dpsi*Spsi*Cthet*Wx - dthet*Cpsi*Sthet*Wx + Cpsi*Cthet*Wxz*dz );
dv = dv - ( ( -dpsi*Spsi*Sthet*Sphi + dthet*Cpsi*Cthet*Sphi + dphi*Cpsi*Sthet*Cphi - dpsi*Cpsi*Cphi + dphi*Spsi*Sphi )*Wx +...
     ( Cpsi*Sthet*Sphi - Spsi*Cphi )*Wxz*dz );
dw = dw - ( ( -dpsi*Spsi*Sthet*Cphi + dthet*Cpsi*Cthet*Cphi - dphi*Cpsi*Sthet*Sphi + dpsi*Cpsi*Sphi + dphi*Spsi*Cphi )*Wx +...
      (Cpsi*Sthet*Cphi + Spsi*Sphi)*Wxz*dz );
  
VT = sqrt(u^2 + v^2 + w^2);
aoa = atan(w/u);
bet = asin(v/VT);
Q = 0.5*1.225*VT*VT;

Ca = cos(aoa); Sa = sin(aoa);
% Cb = cos(bet); Sb = sin(bet);
% Cb = 1; Sb = 0; % when small beta constraint is imposed

Fx = obj.m*(du + w*q-v*r + 9.806*Sthet + Wxz*dz*Cpsi*Cthet);
Fy = obj.m*(dv + r*u-p*w - 9.806*Cthet*Sphi + Wxz*dz*(Cpsi*Sthet*Sphi-Spsi*Cphi));
Fz = obj.m*(dw + p*v-q*u - 9.806*Cthet*Cphi + Wxz*dz*(Cpsi*Sthet*Cphi+Spsi*Sphi));

Mx = (obj.Ix*dpdt-obj.Ixz*drdt + q*r*(obj.Iz-obj.Iy)-p*q*obj.Ixz);
My = (obj.Iy*dqdt + p*r*(obj.Ix-obj.Iz)+(p^2 - r^2)*obj.Ixz);
Mz = (obj.Iz*drdt-obj.Ixz*dpdt + p*q*(obj.Iy-obj.Ix) + q*r*obj.Ixz);

% determination of inputs

a0 = obj.CD0; a1 = obj.CD1; a2 = obj.CD2;

f1 = (Fx*Sa-Fz*Ca)/(Q*obj.S);
f2 = (Fx*Ca+Fz*Sa)/(Q*obj.S);

if abs(Sa)>1e-4
    term1 = (a1*Sa+2*a2*f1*Sa+Ca)/(2*a2*Sa*Sa);
    term2 = (a0+a1*f1+a2*f1*f1+f2)/(a2*Sa*Sa);
    CTx = term1 - sqrt(term1^2-term2);
else
    CTx = (a0+a1*f1+a2*f1*f1+f2);
end

if ~isreal(CTx)
end
    
a0My = obj.Cm0 + obj.Cmalf*aoa + obj.Cmq*(0.5*q*obj.c/VT);
a1My = [obj.Cmdf, obj.Cmde];
a0L =  obj.CL0 + obj.CLalf*aoa + obj.CLq*(0.5*q*obj.c/VT);
a1L =  [obj.CLdf, obj.CLde];

long_CS = [a1My;a1L]\[My/(Q*obj.S*obj.c)-a0My;f1-Sa*CTx-a0L];
df = long_CS(1); de = long_CS(2);

a0l = obj.Clb*bet + obj.Clp*(0.5*p*obj.b/VT) + obj.Clr*(0.5*r*obj.b/VT);
a1l = [obj.Clda, obj.Cldr, 0];

a0n = obj.Cnb*bet + obj.Cnp*(0.5*p*obj.b/VT) + obj.Cnr*(0.5*r*obj.b/VT);
a1n = [obj.Cnda, obj.Cndr, 0];

a0Fy = obj.CYb*bet + obj.CYp*(0.5*p*obj.b/VT) + obj.CYr*(0.5*r*obj.b/VT);
a1Fy = [obj.CYda, obj.CYdr, 0];

a1CTy = [0,0,obj.d/obj.b];

Amat = [Ca*a1l-Sa*a1n;Sa*a1l+Ca*a1n+a1CTy;a1Fy];
bmat = [Mx/(Q*obj.S*obj.b)-(Ca*a0l-Sa*a0n);Mz/(Q*obj.S*obj.b)-(Sa*a0l+Ca*a0n);Fy/(Q*obj.S)-a0Fy];
lat_CS = Amat\bmat;
da = lat_CS(1); dr = lat_CS(2); CTy = lat_CS(3);

Z = [u,v,w,p,q,r,df,da,de,dr,CTx,CTy];
end