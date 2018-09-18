function [c, ceq] = constFun(x, J)
    N = size(J); N = N(1);
    eig_val = x(1:2*N,1);
    eigvec_comp = x(2*N+1:end, 1);
    eig_vec = gen_vec(eigvec_comp, N)
    for i = 1:N
        %ceq(i) = det(J - eig_val(i)*eye(N));
        ceq(i,1) = norm(eig_vec(:,i)) - 1;
    end        
    for i = 1:2:2*N
        ceq(N+1+i:N+1+i+1,1) = (J - complex(eig_val(i,1), eig_val(i+1,1))*eye(N))*eig_vec(:,mod(i,2))
    end
    c = [];
end