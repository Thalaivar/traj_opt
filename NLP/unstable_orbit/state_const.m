function [c, ceq] = state_const(x, p, N, D)
    rk = x(1:N, 1); thetak = x(N+1:2*N,1); T = x(end,1);

    diffmat = [D, zeros(N, N); zeros(N, N), D];
    xdot_cap = -(2/T)*diffmat*[rk;thetak];
    
    fr = zeros(size(rk)); ftheta = zeros(size(thetak));
    for i = 1:length(rk)
        temp = sys_model(0, [rk(i);thetak(i)], p);
        fr(i, 1) = temp(1,1); ftheta(i, 1) = temp(2,1);
    end
    xdot = [fr;ftheta];
    
    c = [];
    ceq = xdot_cap - xdot;
    
    ceq(end+1) = x(1,1) - x(N,1);
    ceq(end+1) = x(N+1,1) - x(2*N,1) + 2*pi;
end