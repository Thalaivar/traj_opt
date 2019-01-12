% aircraft needs N and params set
function f = objfun(X, solution)
    global eigval;
    N = solution.N; n_coeffs = 2*N+1;

    % to be used during trajectory optimisation
    if strcmp(solution.objfun_type, 'traj')
        VR = X(3*n_coeffs+1,1);
        f = VR;
    % to be used when trying new method of stability optimisation
    elseif strcmp(solution.objfun_type, 'floq_new')
        coeffs_x = X(1:n_coeffs,1);
        coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
        coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
        coeffs = [coeffs_x,coeffs_y, coeffs_z];
        tf = X(3*n_coeffs+2,1); VR = X(3*n_coeffs+1,1);
        solution.tf = tf; solution.coeffs = coeffs;
        solution.VR = VR; 

        FTM_expo = get_FTM(solution, 'expo');
        D = eig(FTM_expo); f = 0; param1 = 1; param2 = 10;
        for i = 1:3
             f = f + atan(param1*(abs(D(i)) - 1));
         end
        f = f + param2*VR;
        max_eig = max(abs(D));

        if(isempty(eigval)), eigval = max_eig; end
        if(eigval > max_eig), eigval = max_eig; solution.X = X; end

    end
        
%     % to be used when estimating floquet expo via collocation (probably
%     % obsolete)
%     elseif nargin > 3 && strcmp(type, 'floq_old')
%         n = 3; % dimension of 3DOF system
%         eigval = complex(X(2*n*M+1,1), X(2*n*M+2,1));
%         f = -real(eigval);     
%     end
end