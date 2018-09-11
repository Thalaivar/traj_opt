function ydot = sys_model(t, y, p)
    a = p(1); b = p(2); c = p(3); d = p(4); nu = p(5); omega = p(6);
    x1dot = d*nu*y(1,1) + a*y(1,1)^3;
    x2dot = omega + c*nu + b*y(1,1)^2;
    ydot = [x1dot;x2dot];
end