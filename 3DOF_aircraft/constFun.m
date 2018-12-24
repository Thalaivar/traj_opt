function [c, ceq, dc, dceq] = constFun(X, aircraft, N, M, floq)
    % if eigenvalues are not being estimated
    if(floq == false)
       n_coeffs = 2*N+1;
       coeffs_x = X(1:n_coeffs,1);
       coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
       coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
       coeffs = [coeffs_x,coeffs_y, coeffs_z];
       VR = X(3*n_coeffs+1,1); tf = X(3*n_coeffs+2,1);
       
       M = 50;
       
       c = zeros(12*(M+1),1);
       t = linspace(0, tf, M);
       for i = 1:M
           j = (i-1)*12 + 1;
           % get flat outputs at time t
           sigma = aircraft.get_traj(t(i), tf, coeffs, N);
           z = sigma(3);
           % get states at time t`
           aircraft = aircraft.get_xu(sigma, VR);
           % constrain states and controls
           V = aircraft.x(1); Cl = aircraft.u(1); nu = aircraft.u(2);
           CT = aircraft.u(3);
           % Cl constraints
           c(j,1) = aircraft.limits(1) - Cl;
           c(j+1,1) = Cl - aircraft.limits(2);
           % V contraints
           c(j+2,1) = aircraft.limits(3) - V;
           c(j+3,1) = V - aircraft.limits(4);
           % nu constraints
           c(j+4,1) = aircraft.limits(5) - nu;
           c(j+5,1) = nu - aircraft.limits(6);
           % CT constraints
           c(j+6,1) = aircraft.limits(7) - CT;
           c(j+7,1) = CT - aircraft.limits(8);
           % hmin constraints
           c(j+8,1)  = z + 0.5*aircraft.b*sin(nu);
           c(j+9,1) = z - 0.5*aircraft.b*sin(nu);
           %c(j+10,1) = -z - 100;
           % gamma constraints
           c(j+10,1) = aircraft.limits(9) - aircraft.x(2);
           c(j+11,1) = aircraft.x(2) - aircraft.limits(10);
           
       end
       
       ceq = [];
    
    % if eigenvalues are being estimated
    else      
       coeffs = aircraft.traj_params.coeffs;
       VR = aircraft.traj_params.VR; tf = aircraft.traj_params.tf;
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
           beta(:,i) = eigvec_comp(i,4:6)';
       end
       
       % actual_dot_comp :
       %    size = Mx6
       %    same structure as eigvec_comp but is made of the actual
       %    derivative
       actual_dot_comp = zeros(M,6);
       for i = 1:M
           t = cheb_t(i);
           A = aircraft.get_jac(t, tf, VR, coeffs, N);
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
    if nargout > 2 % gradient of the constraints
      dc = [];
      dceq = [];
    end        
end

