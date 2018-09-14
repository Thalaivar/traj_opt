function [c, ceq] = constFun(x, J)
    N = size(J); N = N(1);
    eig_val = x(1:N,1); eigvec_comp = x(N+1:end, 1);
    eig_vec = gen_vec(eigvec_comp, N);
    for i = 1:N
        %ceq(i) = det(J - eig_val(i)*eye(N));
        ceq(i) = norm(eig_vec(:,i)) - 1;
    end        
    c = [];
end