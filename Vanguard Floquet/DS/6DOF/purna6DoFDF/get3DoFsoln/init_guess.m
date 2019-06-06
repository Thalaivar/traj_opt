function Zinit = init_guess(ac)

% inclined circle guess
N = 50;

tfin = 10; 
t = ((tfin/N)*(1:N))'; % must be a column vector

VR = 0.1;
r = 20;
x = r*(-1+cos(2*pi*t/tfin));       dx = -r*(2*pi/tfin)*sin(2*pi*t/tfin);         ddx = -r*((2*pi/tfin)^2)*cos(2*pi*t/tfin);       
y = -r*sqrt(2)*sin(2*pi*t/tfin);   dy = -sqrt(2)*r*(2*pi/tfin)*cos(2*pi*t/tfin); ddy = sqrt(2)*r*((2*pi/tfin)^2)*sin(2*pi*t/tfin);      
z = -r*(1-cos(2*pi*t/tfin))-0.11;  dz = -r*(2*pi/tfin)*sin(2*pi*t/tfin);         ddz = -r*((2*pi/tfin)^2)*cos(2*pi*t/tfin);  

V(N,1) = 0; chi1(N,1) = 0; gam(N,1) = 0;
CL(N,1) = 0; mu(N,1) = 0;
for j = 1:N

    Wx = VR*(-z(j)).^ac.p_exp;
    Wxz = (ac.p_exp*VR)*((-z(j)).^ac.p_exp)./z(j);
    wind = [Wx,Wxz];

    Z = evalZ_3DOF([x(j),y(j),z(j)],[dx(j),dy(j),dz(j);ddx(j),ddy(j),ddz(j)],wind,ac);
    V(j) = Z(1);
    chi1(j) = Z(2);
    gam(j) = Z(3);
    CL(j) = Z(4);
    mu(j) = Z(5);

end
chi1 = unwrap(chi1);
chi = chi1-(t/(t(end)-t(1)))*(chi1(end)-chi1(1));
chi0 = (chi1(end)-chi1(1))/(t(end)-t(1));
Zinit = [V;chi;gam;x;y;z;CL;mu;chi0;VR;tfin];

% Zinit = expander(Zinit,50);

end