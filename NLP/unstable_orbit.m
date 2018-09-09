N = 100;

y0 = [0.1;0.1];
tspan = [0, 2];
[D, cheb_x] = cheb_diff(N);
cheb_t = ((tspan(2)-tspan(1))/2)*cheb_x + (sum(tspan))/2;

[t, ode_samples] = ode45(@system, cheb_t, y0);

r = y(:,1); theta = y(:,2);

%x = zeros(size(r)); y = zeros(size(r));
%for i = 1:length(r)
%   x(i,1) = r(i, 1)*cos(theta(i, 1)); y(i,1) = r(i, 1)*sin(theta(i, 1));
%end





function ydot = system(t, y)
    r = y(1,1); theta = y(2,1);
    d = -2.0; nu = -pi; a = -3.0; c = pi; b = 5.0; omega = pi;
    
    rdot = d*nu*r + a*(r^3);
    thetadot = omega + c*nu + b*(r^2);
    ydot = [rdot; thetadot];
end