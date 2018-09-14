function J = get_jac(x, p)
    y1 = x(1,1); y2 = x(2,1);
    
    J = [p(1) - p(2)*y2, -p(2)*y1; p(3)*y2, p(3)*y1 - p(4)];
end