function D = analyze_FTM(ac, div_tol)
    %% time evolution
    options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11, 'Events', @(t,y) z_event(t,y));
    data = simulate_traj(ac, 200, options);
    
    %% calculating average values for V and gamma
    t = linspace(0, ac.tf, 1000);
    V_avg = 0;
    for i = 1:length(t)
        sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
        ac = ac.get_xu(sig);
        V_avg = V_avg + ac.x(1);
    end
    V_avg = V_avg/length(t);
    
    %% find point of loss of nonlinearity
    t = data.t; y = data.y;
    n_max = floor(t(end)/ac.tf);
    for i = 1:n_max-1
        [~,j] = min(abs(t - i*ac.tf));
        [~,k] = min(abs(t - (i+1)*ac.tf));
        V_div = abs((y(j,1) - y(k,1))/V_avg);
        gam_div = abs(y(j,3) - y(k,3));
        if(V_div >= div_tol(1) || gam_div >= div_tol(2)),  break; end
    end
    
%% FTM by expo
%     FTM = eye(3);
%     for l = k:-1:j
%         del_t = t(l) - t(l+1);
%         Jk = get_A(ac, t(l), y(l,:)'); Jk = Jk(1:3,1:3);
%         FTM = FTM*expm(Jk*del_t);
%     end
    
%% FTM by time marching
    % state at almost periodic time period
    div_state = y(j:k,:);
    
    % adaptive time march
    y_0 = eye(6); tspan = t(j:k); y_T = zeros(6,6);
    for i = 1:6
        [~,y] = ode45(@(t,y) model(t, y, ac, div_state, tspan), tspan, y_0(:,i));
        y_T(:,i) = y(end,:)';
    end
    FTM = y_T(1:3, 1:3);
    
     % fixed step time march (RK4)
%     y_0 = eye(6); tspan = t(j:k); y_T = zeros(6,6);
%     for i = 1:6
%         [~,y] = RK4(@(t,y) model(t, y, ac, div_state, tspan), [tspan(1), tspan(end)], y_0(:,i), 1e-3);
%         y_T(:,i) = y(end,:)';
%     end
%     FTM = y_T(1:3, 1:3);

    function ydot = model(t, y, ac, div_state, tspan)
        state = interp1(tspan, div_state, t);
        X = state';
        A = get_A(ac, t, X);
            ydot = A*y;
    end    

    D = eig(FTM);
end