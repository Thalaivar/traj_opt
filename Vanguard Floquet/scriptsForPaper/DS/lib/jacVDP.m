function A = jacVDP(X, trajData)
    mu = trajData.mu;
%     N = 0.5*(length(X)-1);
    A = zeros(2);
    x = X(1); xdot = X(2);
    A(1,2) = 1;
    A(2,1) = -1 - 2*mu*x*xdot; A(2,2) = mu*(1-x^2);
end