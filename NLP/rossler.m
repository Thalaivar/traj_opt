N = 20;

% simulate Rossler system
y0 = [1;2;3];
tspan = [0, 100];
[t, y] = ode45(@RS, tspan, y0);
plot(t, y(:,1), '-r', t, y(:,2), '-b', t, y(:,3), '-m');

% to simulate the system
function ydot = RS(t, y)
    x1 = y(1,1); x2 = y(2,1); x3 = y(3,1);
    p = [0.2, 0.2, 5.7];
    x1dot = -1*x2 - x3;
    x2dot = x1 + p(1)*x2;
    x3dot = p(2) + x3*(x1 - p(3));
    ydot = [x1dot;x2dot;x3dot];
end