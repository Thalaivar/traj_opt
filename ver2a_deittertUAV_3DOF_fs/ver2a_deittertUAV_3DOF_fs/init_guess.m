function [Xinit,N,M] = init_guess()

M = 8; % no. of harmonics
N = 50; % no. of points where cineq and ceq are imposed

tf = 10; % 20 for "00"
t = linspace(0,tf,3*N)';
VR = 0.1;
r = 20;

% "O"
x = r*(-1+cos(2*pi*t/tf));
y = -r*sqrt(2)*sin(2*pi*t/tf);  
z = -r*(1-cos(2*pi*t/tf))-0.11;

% % "00" needs a larger tf (typical: tr = 20 and VR = 0.1)
% x = -r*sin(4*pi*t/tf);                  
% y = -2*r*sqrt(2)*sin(2*pi*t/tf);        
% z = -r*(sin(4*pi*t/tf)+1)-0.11;  

 

fourCoeff = fitFS([x,y,z],tf,M);

% load pseudoSpec_loiter_001_fourCoeff
% tf = 11.9998;
% VR = 0.0931;
% M = 8;


Xinit(3*(2*M+1)+2,1) = 0;
Xinit(end-1) = tf;
Xinit(end) = VR;
for i = 1:3
    Xinit((i-1)*(2*M+1)+1:i*(2*M+1)) = fourCoeff(i,:);
end



end