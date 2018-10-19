function [c, ceq] = constFun_floq(X, limits, model_par, N, M)
    % limits is of the form:
    %       limits = [Clmax, Vmax, nu_min, nu_max, Tmin, Tmax, hmin]    
    
    % model_par is of the form:
    %       model_par = [m, rho, S, g,Cd0, Cd1, Cd2, b]
    
    % z is of the form:
    %       z = [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot]
    
    % n1 : n_coeffs; n2 : n_phase_angles; n3 : n_eigvec_comp;
    n1 = N+1; n2 = N; n3 = 6*M;
    a_h = X(1:n1,1); a_x = X(n1+1:2*n1,1); a_y = X(2*n1 + 1:3*n1,1);
    eta_h = X(3*n1+1:3*n1+n2,1); eta_x = X((3*n1+n2+1):(3*n1+2*n2),1); eta_y = X((3*n1+2*n2+1):3*(n1+n2),1);
    ph = [a_h,[eta_h;0]]; px = [a_x,[eta_x;0]]; py = [a_y,[eta_y;0]]; VR = X(3*(n1+n2)+1,1); tf = X(3*(n1+n2)+2,1);
    
    wind_par = [VR, 1];
    
    % construct eigenvector at diff time points
    eig_vec_comp1 = X(3*(n1+n2)+2+1:3*(n1+n2)+2+n3,1); eigval = X(end,1);
    eig_vec = zeros(6,M);
    for i = 1:M
        j = (i-1)*6 + 1;
        eig_vec(:,i) = eig_vec_comp1(j:j+5,1);
    end 
    % matrix with columns as the values of the components of the
    % eigenvectors at diff time points
    eig_vec_comp2 = eig_vec'; 
    
    [D, cheb_x] = cheb_diff(M-1);
    cheb_t = 0.5*tf*(1-cheb_x);
    
    % cheb derivative
    zero_mat = zeros(M);
    diff_mat = [D       , zero_mat, zero_mat, zero_mat, zero_mat, zero_mat;
                zero_mat, D       , zero_mat, zero_mat, zero_mat, zero_mat;
                zero_mat, zero_mat, D       , zero_mat, zero_mat, zero_mat;
                zero_mat, zero_mat, zero_mat, D       , zero_mat, zero_mat;
                zero_mat, zero_mat, zero_mat, zero_mat, D       , zero_mat;
                zero_mat, zero_mat, zero_mat, zero_mat, zero_mat, D       ;];
    
    u1 = eig_vec_comp2(:,1); u5 = eig_vec_comp2(:,5);
    u2 = eig_vec_comp2(:,2); u4 = eig_vec_comp2(:,4); 
    u3 = eig_vec_comp2(:,3); u6 = eig_vec_comp2(:,6);
    udot_cap = -(2/tf)*diff_mat*[u1;u2;u3;u4;u5;u6];
    
    % actual derivative
    udot = zeros(6*M,1);
    a = [a_h, a_x, a_y]; eta = [eta_h, eta_x, eta_y];
    for i = 1:M
        temp_A = get_jac(cheb_t(i), a, eta, tf, model_par, wind_par); 
        temp_udot = (temp_A - eigval*eye(6))*eig_vec(:,i);
        for k = 1:6
            udot((k-1)*M + i,1) = temp_udot(k,1);
        end
    end
    
    % constraints for derivative and norm
    ceq = zeros(6*M + M);
    % diff constraints
    ceq(1:6*M,1) = udot - udot_cap;
    
    % magnitutde of eigenvectors is 1
    for i = 1:M
        ceq(6*M+i) = norm(eig_vec(:,i)) - 1;
    end
    
    O = 6*N;
    c = zeros(8*O,1);
    for i = 1:O+1
        j = (i-1)*8 + 1;
        t = (i-1)*tf/O;
        [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot] = get_traj(t, ph, px, py, tf);
        z = [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot];
        c(j:j+7,1) = ineq_constr(z, limits, wind_par, model_par);
    end
end