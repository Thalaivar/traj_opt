function new_guess = mod_init_guess(aircraft, X, M, N)
   eigval = complex(X(end-1), X(end));
   eigvec_comp = [X(1:6*M,1), X(6*M+1:2*6*M,1)];
   alpha_prev = zeros(6,M); beta_prev = zeros(6,M);
   for i = 1:M
       j = (i-1)*6 + 1;
       alpha_prev(:,i) = eigvec_comp(j:j+5,1);
       beta_prev(:,i)  = eigvec_comp(j:j+5,2);
   end
   
   tf = aircraft.traj_params.tf; VR = aircraft.traj_params.VR; coeffs = aircraft.traj_params.coeffs;
   [~,cheb_x] = cheb_diff(M-1); t_prev = 0.5*tf*(1-cheb_x);
   % add points between the two original points to generate a finer mesh
   alpha_new = zeros(6,(M-1)); beta_new = zeros(6,(M-1)); 
   for i = 1:M-1
       % time point considered is midpoint of eexisting two points
       t = (t_prev(i) + t_prev(i+1))/2; hk = t_prev(i+1) - t_prev(i);
       % jacobian at new point
       A = aircraft.get_jac(t, tf, VR, coeffs, N); I = eye(6);
       % hermite-simpson interpolation
       ak1 = alpha_prev(:,i); ak2 = alpha_prev(:,i+1);
       bk1 = beta_prev(:,i); bk2 = beta_prev(:,i+1);
       fka_1 = (A - real(eigval)*I)*ak1 + imag(eigval)*bk1;
       fka_2 = (A - real(eigval)*I)*ak2 + imag(eigval)*bk2;
       fkb_1 = (A - real(eigval)*I)*bk1 - imag(eigval)*ak1; 
       fkb_2 = (A - real(eigval)*I)*bk2 - imag(eigval)*ak2;
       alpha_new(:,i) = 0.5*(ak1 + ak2) + (hk/8)*(fka_1 - fka_2);
       beta_new(:,i)  = 0.5*(bk1 + bk2) + (hk/8)*(fkb_1 - fkb_2);
   end
   
   new_guess = zeros(12*M+12*(M-1),1);
   for i = 1:M-1
       j = (i-1)*12 + 1;
       new_guess(j:j+11,1) = [alpha_prev(:,i);alpha_new(:,i);];
       new_guess(6*M+6*(M-1)+j:6*M+6*(M-1)+j+11,1) = [beta_prev(:,i);beta_new(:,i);];
   end
   new_guess(12*(M-1)+1:6*(M-1)+6*M,1) = alpha_prev(:,M);
   new_guess(6*M+6*(M-1)+12*(M-1)+1:12*M+12*(M-1),1) = beta_prev(:,M);
   new_guess(end+1,1) = real(eigval); new_guess(end+1,1) = imag(eigval);
end