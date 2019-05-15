function Zinit = init_guess(ac)

% load('../3DoF_DS_dfFourier_linear/results/eightLin_001');
% Zinit = [dsSoln.V;
%          dsSoln.chi;
%          dsSoln.gam;
%          dsSoln.x;
%          dsSoln.y;
%          dsSoln.z;
%          dsSoln.CL;
%          dsSoln.mu;
%          dsSoln.VR;
%          dsSoln.T];
% N = dsSoln.N;

% inclined double-circle guess
N = 50;

tfin = 20; 
t = ((tfin/N)*(1:N))'; % must be a column vector

VR = 0.09;
r = 20;
x = -r*sin(4*pi*t/tfin);           dx = -(4*pi*r/tfin)*cos(4*pi*t/tfin);         ddx = ((4*pi/tfin)^2)*r*sin(4*pi*t/tfin);       
y = -2*r*sqrt(2)*sin(2*pi*t/tfin); 
dy = -4*r*sqrt(2)*(pi/tfin)*cos(2*pi*t/tfin); 
ddy = 8*r*sqrt(2)*((pi/tfin)^2)*sin(2*pi*t/tfin);      
z = -r*(sin(4*pi*t/tfin)+1)-0.11;  dz = -(4*pi*r/tfin)*cos(4*pi*t/tfin);         ddz = ((4*pi/tfin)^2)*r*sin(4*pi*t/tfin);  

V(N,1) = 0; chi(N,1) = 0; gam(N,1) = 0;
CL(N,1) = 0; mu(N,1) = 0;
for j = 1:N

    Wx = VR*(-z(j)).^ac.p_exp;
    Wxz = (ac.p_exp*VR)*((-z(j)).^ac.p_exp)./z(j);
    wind = [Wx,Wxz];

    Z = evalZ_3DOF([x(j),y(j),z(j)],[dx(j),dy(j),dz(j);ddx(j),ddy(j),ddz(j)],wind,ac);
    V(j) = Z(1);
    chi(j) = Z(2);
    gam(j) = Z(3);
    CL(j) = Z(4);
    mu(j) = Z(5);

end
Zinit = [V;chi;gam;x;y;z;CL;mu;VR;tfin];

% Zinit = expander(Zinit,50);

end