function A = jacGlycolytic(X, trajData)
    a = trajData.a;
    x = X(1); y = X(2);
    A(1,1) = -1 + 2*x*y;
    A(1,2) = a + x^2;
    A(2,1) = -2*x*y;
    A(2,2) = -a - x^2;
end