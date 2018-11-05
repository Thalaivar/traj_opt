
function [c, ceq, dc, dceq] = constFun(X, aircraft, N, M, floq)
    % if eigenvalues are not being estimated
    if(floq == false)
       n_coeffs = 2*N+1;
       coeffs_x = X(1:n_coeffs,1);
       coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
       coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
       coeffs = [coeffs_x,coeffs_y, coeffs_z];
       VR = X(3*n_coeffs+1,1); tf = X(3*n_coeffs+2,1);
    
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
       eigvec_comp = [X(1:3*M,1), X(3*M+1:2*3*M,1)]; eigval = complex(X(end-1,1), X(end,1));
       
       c = [];

       % alpha : vector made up of real part of eigvec
       % beta  : vector made up of img part of eigvec
       % generate time series of components of eigenvectors at diff time points.
       alpha = zeros(3,M); beta = zeros(3,M);
       for i = 1:M
           j = (i-1)*3 + 1;
           alpha(:,i) = eigvec_comp(j:j+2,1);
           beta(:,i)  = eigvec_comp(j:j+2,2);
       end
       % matrix with colums as time series of each component
       eigvec_comp_real = alpha'; eigvec_comp_img = beta';
        
       ceq = zeros(2*3*M,1);
       % cheb diff matrix for 6 components of eigvector
       [D, cheb_x] = cheb_diff(M-1);
       diffmat = [D       ,zeros(M),zeros(M);
                  zeros(M),D       ,zeros(M);
                  zeros(M),zeros(M),D       ];
              
       % vector of the form [u1', u2', .. , u6']', where ui is 
       % a column vector of  the value of i'th component of the eigenvector
       % at M points
       udot_real = zeros(3*M,1); udot_img = zeros(3*M,1);
       for i = 1:3
           j = (i-1)*M+1;
           udot_real(j:j+M-1,1) = eigvec_comp_real(:,i);
           udot_img(j:j+M-1,1) = eigvec_comp_img(:,i);
       end
       % cheb interpolated derivative
       alphadot_cap = -(2/tf)*diffmat*udot_real;
       betadot_cap  = -(2/tf)*diffmat*udot_img;
       
       % calculating actual derivative
       cheb_t = 0.5*tf*(1-cheb_x);
       % matrix has columns as time series of components of eigenvector at
       % M points (actual value)
       udot_real_actual = zeros(M,3); udot_img_actual = zeros(M,3);
       for i = 1:M
           I = eye(3);
           A = aircraft.get_jac(cheb_t(i), tf, VR, coeffs, N);
           temp_udot_real = (A - real(eigval)*I)*alpha(:,i) + imag(eigval)*beta(:,i);
           temp_udot_img  = (A - real(eigval)*I)*beta(:,i)  - imag(eigval)*alpha(:,i); 
           udot_real_actual(i,:) = temp_udot_real'; 
           udot_img_actual(i,:) = temp_udot_img';
       end
       alphadot = zeros(3*M,1); betadot = zeros(3*M,1);
       for i = 1:3
           j = (i-1)*M+1;
           alphadot(j:j+M-1,1) = udot_real_actual(:,i);
           betadot(j:j+M-1,1) = udot_img_actual(:,i);
       end
       ceq(1:3*M,1) = alphadot_cap - alphadot;
       ceq(3*M+1:6*M,1) = betadot_cap - betadot;
       
%        % norm of eigenvectors constrained to be = 1
%        eigvec = complex(alpha, beta);
%        for i = 1:M
%            ceq(12*M+i,1) = norm(eigvec(:,i)) - 1;
%        end
       
    end
    if nargout > 2 % gradient of the constraints
      dc = [];
      dceq = [];
    end        
end

