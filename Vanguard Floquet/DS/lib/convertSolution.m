function X = convertSolution(p, N, sol, type)
    if strcmp(type, 'coeffs')
        n_h = (((length(sol)-2)-3)/3)/2;
        n_c = 1+n_h*2;
        T = sol(n_c*3+2); VR = sol(n_c*3+1);
        coeffs = [sol(1:n_c), sol(n_c+1:2*n_c), sol(2*n_c+1:3*n_c)];
        [~,t] = fourierdiff(N);
        t = T*t/(2*pi);
        x = zeros(N,6); u = zeros(N,3);
        for i = 1:N
            [xx,uu] = get_traj(t(i), T, coeffs, n_h, VR, p);
            x(i,:) = xx; u(i,:) = uu;
        end
        x(:,2) = unwrap(x(:,2));
        X = zeros(9*N+2,1);
        for i = 1:9
            j = (i-1)*N;
            if(i <= 6)
                X(j+1:j+N) = x(:,i);
            else
                X(j+1:j+N) = u(:,i-6);
            end
        end
        X(9*N+1) = T; X(9*N+2) = VR; 
    elseif strcmp(type, 'xyz')
        NN = (length(sol) - 2)/3;
        xyz = [sol(1:NN), sol(NN+1:2*NN), sol(2*NN+1:3*NN)];
        T = sol(3*NN+1); VR = sol(3*NN+2);
        coeffs = fourierbasis(xyz, T, 8);
        df_vec = [coeffs(:,1); coeffs(:,2); coeffs(:,3); VR; T];
        X = df_to_fourier_colloc(p, N, df_vec, 'coeffs');
    end
end