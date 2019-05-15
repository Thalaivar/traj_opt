function Z = evalZ_3DOF(Y,dY,wind,obj)
% 27/12
% Z: non-flat states and controls
% evaluates Z = [V,chi,gam,CL,mu,CT]
% the flat outputs and its first two derivatives.
% Y = [x,y,z]
% Y: 1x3
% dY: 2x3

% x = Y(1); y = Y(2); z = Y(3);

dx = dY(1,1); dy = dY(1,2); dz = dY(1,3);
ddx = dY(2,1); ddy = dY(2,2); ddz = dY(2,3);

Wx = wind(1); Wxz = wind(2);

V = sqrt( (dx-Wx)^2 + dy^2 + dz^2 );

chi = atan2(dy,(dx-Wx));
% if (chi-init_heading)<0
%     chi = chi+2*pi;
% end    

gam = -asin(dz/V);

dV = (dx*ddx - dx*dz*Wxz - ddx*Wx + Wx*Wxz*dz + dy*ddy + dz*ddz)/V;
dchi = (dx*ddy - ddy*Wx - dy*ddx + dy*dz*Wxz)/(dy^2 + dx^2 + Wx^2 - 2*dx*Wx);
dgam = (dV*dz - V*ddz)/(V*sqrt(V^2 - dz^2));

Cchi = cos(chi); Schi = sin(chi);
Cgam = cos(gam); Sgam = sin(gam);

Q = 0.5*1.225*V*V;

mu = atan( ( V*Cgam*dchi - Wxz*dz*Schi )/( V*dgam + 9.806*Cgam - Wxz*Cchi*Sgam*dz ) );
CL = (obj.m*V*Cgam*dchi - obj.m*Wxz*dz*Schi)/(Q*obj.S*sin(mu));

D = Q*obj.S*(obj.CD0 + obj.CD1*CL + obj.CD2*CL*CL);

T = obj.m*dV + D + obj.m*9.806*Sgam + obj.m*Wxz*dz*Cchi*Cgam;
CT = (T/(Q*obj.S));

Z = [V,chi,gam,CL,mu,CT];
% dX = [dV,dchi,dgam];

end
