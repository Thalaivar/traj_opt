function xdot = vdp(x)
    nu = 0.2;
    x1 = x(1,1); x2 = x(2,1);
    
    x1dot = x2;
    x2dot = nu*(1 - x1^2)*x2 - x1;
    
    xdot = [x1dot; x2dot];
end