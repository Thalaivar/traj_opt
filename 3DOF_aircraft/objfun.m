% aircraft needs N and params set
function f = objfun(X, aircraft, type, M)
    N = aircraft.N; n_coeffs = 2*N+1;
    if nargin <= 3
        % to be used during trajectory optimisation
        if strcmp(type, 'traj')
            VR = X(3*n_coeffs+1,1);
            f = VR;
        % to be used when trying new method of stability optimisation
        elseif strcmp(type, 'floq_new')
            coeffs_x = X(1:n_coeffs,1);
            coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
            coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
            coeffs = [coeffs_x,coeffs_y, coeffs_z];
            VR = X(3*n_coeffs+1,1); tf = X(3*n_coeffs+2,1);

            aircraft.tf = tf; aircraft.VR = VR;
            aircraft.coeffs = coeffs;

            [FTM_expo, ~] = get_FTM(aircraft);
            D = (1/tf)*log(eig(FTM_expo));
            f = -max(real(D));
        end
    % to be used when estimating floquet expo via collocation (probably
    % obsolete)
    elseif nargin > 3 && strcmp(type, 'floq_old')
        n = 3; % dimension of 3DOF system
        eigval = complex(X(2*n*M+1,1), X(2*n*M+2,1));
        f = -real(eigval);
    end
end