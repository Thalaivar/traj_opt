function ydot = sys_model(t, y, p)
    y1 = y(1,1); y2 = y(2,1);
    ydot1 = p(1)*y1 - p(2)*y1*y2;
    ydot2 = p(3)*y1*y2 - p(4)*y2;
    ydot = [ydot1; ydot2];
end