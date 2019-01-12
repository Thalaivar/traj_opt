% to be used when estimating eigenvalues via pseudospectral collocation
% (probably obsolete)
% aircraft object needs to have N, tf, VR, coeffs and the params set
% NOTE: if init_guess is used, M should match
function [c, ceq] = constFun_floq(X, aircraft, M)
       % eigvec_comp : 
       %    size = Mx6
       %    columns : alpha_1, alpha_2, ... , beta_3 (alpha_i is i'th component of alpha)
       %    rows : i'th row is i'th collocation point            
       eigvec_comp = zeros(M,6);
       for i = 1:6
           j = (i-1)*M + 1;
           eigvec_comp(:,i) = X(j:j+M-1,1)';
       end
       
       eigval = complex(X(end-1,1), X(end,1));
       
       c = [];
 
       % cheb diff matrix for 6 components of eigvector
       [D, cheb_x] = cheb_diff(M-1);
              
       % dot_cap : 
       %    size = 6Mx1
       %    is of the form [alpha_1_dot_cap', ... , beta_3_dot_cap']'
       %    alpha_1_dot_cap represents the cheb derivative of alpha_1
       dot_cap = zeros(6*M,1);
       for i = 1:6
            j = (i-1)*M + 1;
            dot_cap(j:j+M-1,1) = (-2/tf)*D*eigvec_comp(:,i);
       end
       
       % calculating actual derivative
       cheb_t = 0.5*tf*(1-cheb_x);
       
       % alpha, beta:
       %    size = 3xM
       %    each column is value of alpha at that collocation point
       alpha = zeros(3,M); beta = zeros(3,M);
       for i = 1:M
           alpha(:,i) = eigvec_comp(i,1:3)';
           beta(:,i)  = eigvec_comp(i,4:6)';
       end
       
       % actual_dot_comp :
       %    size = Mx6
       %    same structure as eigvec_comp but is made of the actual
       %    derivative
       actual_dot_comp = zeros(M,6);
       for i = 1:M
           t = cheb_t(i);
           A = aircraft.get_jac(t);
           alpha_dot = (A - real(eigval)*eye(3))*alpha(:,i) + imag(eigval)*beta(:,i);
           beta_dot  = (A - real(eigval)*eye(3))*beta(:,i)  - imag(eigval)*alpha(:,i);
           actual_dot_comp(i,1:3) = alpha_dot';
           actual_dot_comp(i,4:6) = beta_dot';
       end
       
       % dot : 
       %    size = 6Mx1
       %    is of the same for as dot_cap, except the values are now real derivatives
       dot = zeros(6*M,1);
       for i = 1:6
            j = (i-1)*M + 1;
            dot(j:j+M-1,1) = actual_dot_comp(:,i);
       end
       
       ceq = dot - dot_cap; 
       ceq(end+1:end+6,1) = eigvec_comp(1,:)' - eigvec_comp(end,:)'; 
end

