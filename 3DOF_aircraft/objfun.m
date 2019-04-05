% aircraft needs N and params set
function f = objfun(X, aircraft, type, params, stab_type, M)
    N = aircraft.N; n_coeffs = 2*N+1;
    if nargin <= 5
        % to be used during trajectory optimisation
        if strcmp(type, 'traj')
            VR = X(3*n_coeffs+1,1);
            f = VR;
        
        elseif strcmp(type, 'energy')
            coeffs_x = X(1:n_coeffs,1);
            coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
            coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
            coeffs = [coeffs_x,coeffs_y, coeffs_z];
            tf = X(3*n_coeffs+1,1);
            aircraft.tf = tf; aircraft.coeffs = coeffs;
            
            t = linspace(0, aircraft.tf, 1000);
            state = zeros(length(t), 6);
            E_a = zeros(1, length(t));
            for i = 1:length(t)
                sig = get_traj(t(i), aircraft.tf, aircraft.coeffs, aircraft.N);
                aircraft = aircraft.get_xu(sig);
                state(i,:) = [aircraft.x(1), aircraft.x(3), aircraft.x(2), sig(1), sig(2), sig(3)];
                E_a(i) = 0.5*aircraft.m*(state(i,1)^2) + aircraft.m*aircraft.g*(-state(i,6));
            end
            
            f = -mean(E_a);
            
        % to be used when trying new method of stability optimisation
        elseif strcmp(type, 'stability')
            coeffs_x = X(1:n_coeffs,1);
            coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
            coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
            coeffs = [coeffs_x,coeffs_y, coeffs_z];
            tf = X(3*n_coeffs+1,1);
            aircraft.tf = tf; aircraft.coeffs = coeffs;
           
            FTM_expo = get_FTM(aircraft, 'friedmann');
            f = 0; %param1 = 0.1; param2 = 1; % 10, 0.1
            if(aircraft.p ~= 1)
                FTM_expo = [FTM_expo(1:3,1:3),FTM_expo(1:3,6);FTM_expo(6,1:3),FTM_expo(6,6)];
                D = eig(FTM_expo);
                f = 10*(max(abs(D)) - 1)^2;
                for i = 1:4
                    if strcmp(stab_type, 'stable')
                         f = f + atan(params(1)*(abs(D(i))- 1));
                    elseif strcmp(stab_type, 'unstable')
                        f = f + atan(params(1)*((abs(D(i)))^(-1)- 1));
                    end
                end
            else
                FTM_expo = FTM_expo(1:3,1:3);
                D = eig(FTM_expo);
                for i = 1:3
                     f = f + atan(params(1)*(abs(D(i))- 1));
                    % f = f + atan(param1*((abs(D(i)))^(-1)- 1));
                end
            end
            
        elseif strcmp(type, 'stab_traj')
            coeffs_x = X(1:n_coeffs,1);
            coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
            coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
            coeffs = [coeffs_x,coeffs_y, coeffs_z];
            tf = X(3*n_coeffs+2,1); VR = X(3*n_coeffs+1,1);
            aircraft.tf = tf; aircraft.coeffs = coeffs;
            aircraft.VR = VR;
           
            FTM_expo = get_FTM(aircraft, 'friedmann');
            f = 0; %param1 = 0.1; param2 = 1; % 10, 0.1
            if(aircraft.p ~= 1)
                FTM_expo = [FTM_expo(1:3,1:3),FTM_expo(1:3,6);FTM_expo(6,1:3),FTM_expo(6,6)];
                D = eig(FTM_expo);
                for i = 1:4
                    if strcmp(stab_type, 'stable')
                        f = f + atan(params(1)*(abs(D(i))- 1));
                    elseif strcmp(stab_type, 'unstable')
                        f = f + atan(params(1)*((abs(D(i)))^(-1)- 1));
                    end
                end
            else
                FTM_expo = FTM_expo(1:3,1:3);
                D = eig(FTM_expo);
                for i = 1:3
                     f = f + atan(params(1)*(abs(D(i))- 1));
                    % f = f + atan(param1*((abs(D(i)))^(-1)- 1));
                end
            end
            f = f + params(2)*aircraft.VR;
        
        end
        
    % to be used when estimating floquet expo via collocation (probably
    % obsolete)
    elseif nargin > 5 && strcmp(type, 'floq_old')
        n = 3; % dimension of 3DOF system
        eigval = complex(X(2*n*M+1,1), X(2*n*M+2,1));
        f = -real(eigval);     
    end
end