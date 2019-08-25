function K = friedmann_K(h, t, ac)
    A_psi = ac.get_jac(t, 'FD');
    A_psi_h_by_2 = ac.get_jac(t+0.5*h, 'FD');
    A_psi_h = ac.get_jac(t+h, 'FD');
    E = A_psi_h_by_2*(eye(6) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(6) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(6) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(6) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end