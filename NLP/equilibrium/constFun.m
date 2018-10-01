function [c , ceq] = constFun(x, p)
    N = p(1); eigval_prev = p(2);
    eig_val = x(1,1);
    v = (eig_val - eigval_prev);
    c = -v'*v;
    ceq = [];
end