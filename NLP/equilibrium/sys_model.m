function ydot = sys_model(t, y)
    y1 = y(1,1); y2 = y(2,1);
    y1dot = 1 + (y1^2)*y2 - 4*y1;
    y2dot = 3*y1 - y2*(y1^2);
    ydot = [y1dot;y2dot];
end