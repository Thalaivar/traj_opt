function J = get_jac(x)
    y1 = x(1,1); y2 = x(2,1);
    J = [(2*y1*y2 - 4), y1^2 ; (3 - 2*y1*y2) , -1*y1^2];
end