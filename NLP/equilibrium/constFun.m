function [c, ceq] = constFun(x, J)
    N = size(J); N = N(1);
    eig_val = x(1:N,1); eigvec_comp = x(N+1:end, 1);
    j = 1;
    for i = 1:N
        j = 2*i - 1;
        eig_vec = [eigvec_comp(j); eigvec_comp(j+1)];
        %ceq(i) = det(J - eig_val(i)*eye(N));
        ceq(i) = norm(eig_vec) - 1;
    end        
    c = [];
end