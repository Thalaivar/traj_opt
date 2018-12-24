function f = objfun_floq(X, aircraft, N )
    
    n_coeffs = 2*N+1;
    coeffs_x = X(1:n_coeffs,1);
    coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
    coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
    coeffs = [coeffs_x,coeffs_y, coeffs_z];
    VR = X(3*n_coeffs+1,1); tf = X(3*n_coeffs+2,1);
    
    aircraft.traj_params.tf = tf; aircraft.traj_params.VR = VR;
    aircraft.traj_params.coeffs = coeffs;

    t = linspace(tf, 0, 1000);
    Jk_0 = (aircraft.get_jac(t(1), tf, VR, coeffs, N) - aircraft.get_jac(t(2), tf, VR, coeffs, N))/2;
    H2 = expm(Jk_0*(t(1)-t(2)));
    for i = 2:length(t)-1
        Jk = (aircraft.get_jac(t(i), tf, VR, coeffs, N) + aircraft.get_jac(t(i+1), tf, VR, coeffs, N))/2;
        H2 = H2*expm(Jk*(t(i)-t(i+1)));
    end
    
    D = (1/tf)*log(eig(H2));
    f = -max(real(D));
    
end