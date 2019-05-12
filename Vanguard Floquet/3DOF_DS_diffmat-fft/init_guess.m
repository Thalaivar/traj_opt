function [Xinit,N,which_diff] = init_guess()

% which_diff = 1; % use diffmat
which_diff = 0; % use fft

% N must be even ; no. of points where cineq and ceq are imposed

% % "O"
% % activate all the rate constraints
% % both CT=0 and |CT|<del works
% % interior point is preferred

N = 36;
tf = 10; % 10 for 'O', 20 for '8'
t = ((tf/N)*(1:N))'; % must be a column vector
VR = 0.1; 
r = 20;
x = r*(-1+cos(2*pi*t/tf));
y = -r*sqrt(2)*sin(2*pi*t/tf);  
z = -r*(1-cos(2*pi*t/tf))-0.11;

% % "00" 
% % needs a larger tf (~20)
% % activate only the constraint on the second derivative of mu
% % CT=0 is preferred
% % sqp preferred

% N = 36;
% tf = 20; % 10 for 'O', 20 for '8'
% t = ((tf/N)*(1:N))'; % must be a column vector
% VR = 0.09; 
% r = 20;
% x = -r*sin(4*pi*t/tf);                  
% y = -2*r*sqrt(2)*sin(2*pi*t/tf);        
% z = -r*(sin(4*pi*t/tf)+1)-0.11;  

Xinit(3*N+2,1) = 0;
Xinit(1:end-2) = vertcat(x,y,z);
Xinit(end-1) = tf;
Xinit(end) = VR;


end