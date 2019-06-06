function Zinit = init_guess(ac)

% N = 50;
% tfin = 10;
% t = ((tfin/N)*(1:N))'; fac = 2*pi/tfin; fac2 = fac*fac;
% [D,D2,~] = fourierDiff2(N);
% VR = 0.1;
% r = 20;
% x = r*(-1+cos(2*pi*t/tfin));       dx = fac*D*x;   ddx = fac2*D2*x;       
% y = -r*sqrt(2)*sin(2*pi*t/tfin);   dy = fac*D*y;   ddy = fac2*D2*y;      
% z = -r*(1-cos(2*pi*t/tfin))-0.11;  dz = fac*D*z;   ddz = fac2*D2*z;  
% 
% Psi(N,1) = 0; Thet(N,1) = 0; Phi(N,1) = 0;
% for j = 1:N
% 
%     Wx = VR*(-z(j)).^ac.p_exp;
%     Wxz = (ac.p_exp*VR)*((-z(j)).^ac.p_exp)./z(j);
%     wind = [Wx,Wxz];
% 
%     Z = evalZ_3DOF([x(j),y(j),z(j)],[dx(j),dy(j),dz(j);ddx(j),ddy(j),ddz(j)],wind,ac);
%     Psi(j) = Z(2);
%     Thet(j) = Z(3);
%     Phi(j) = Z(5);
% 
% end
% Psi = unwrap(Psi);
% Psi1 = Psi-(t/(t(end)-t(1)))*(Psi(end)-Psi(1));
% Psi0 = (Psi(end)-Psi(1))/(t(end)-t(1));
% 
% Zinit = [x;y;z;Phi;Thet;Psi1;Psi0;VR;tfin];

% 3DoF DS soln
load('./get3DoFsoln/results/loiterExp_002')
Zinit = [dsSoln.x;dsSoln.y;dsSoln.z;dsSoln.mu;dsSoln.gam;dsSoln.chi;dsSoln.chi0;dsSoln.VR;dsSoln.T];

% load('./get3DoFsoln/results/eightExp_001')
% Zinit = [dsSoln.x;dsSoln.y;dsSoln.z;dsSoln.mu;dsSoln.gam;dsSoln.chi;0;dsSoln.VR;dsSoln.T];

% Zinit = expander(Zinit,50);

end