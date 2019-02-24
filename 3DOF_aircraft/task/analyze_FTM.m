function D = analyze_FTM(ac, div_tol)
    % time evolution
    options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11, 'Events', @(t,y) z_event(t,y));
    data = simulate_traj(ac, 200, options);
    
    % calculating average values for V and gamma
    t = linspace(0, ac.tf, 1000);
    V_avg = 0;
    for i = 1:length(t)
        sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
        ac = ac.get_xu(sig);
        V_avg = V_avg + ac.x(1);
    end
    V_avg = V_avg/length(t);
    
    % find point of loss of nonlinearity
    t = data.t; y = data.y;
    n_max = floor(t(end)/ac.tf);
    max_D = zeros(1,n_max-1);
    for i = 1:n_max-1
        [~,j] = min(abs(t - i*ac.tf));
        [~,k] = min(abs(t - (i+1)*ac.tf));
        FTM = eye(3);
        for l = j:k
            del_t = t(l+1) - t(l);
            Jk = get_A(ac, t(l), y(l,:)'); Jk = Jk(1:3,1:3);
            FTM = FTM*expm(Jk*del_t);
        end
        
        max_D(i) = max(abs(eig(FTM)));
%         V_div = abs((y(j,1) - y(k,1))/V_avg);
%         gam_div = abs(y(j,3) - y(k,3));
%         if(V_div >= div_tol(1) || gam_div >= div_tol(2)),  break; end
    end    
    
    % calulate FTM
    % FTM by exponentials method
%     FTM = eye(3);
%     for i = j:k
%         del_t = t(i+1) - t(i);
%         Jk = get_A(ac, t(i), y(i,:)'); Jk = Jk(1:3,1:3);
%         FTM = FTM*expm(Jk*del_t);
%     end
    
    D = eig(FTM);
end