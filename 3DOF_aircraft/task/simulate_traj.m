function [t, dev_data] = simulate_traj(ac, n, options)
    % time evolution
    tspan = [0, n*ac.tf];
    sig_0 = get_traj(tspan(1), ac.tf, ac.coeffs, ac.N);
    ac = ac.get_xu(sig_0);
    y0 = [ac.x(1);ac.x(3);ac.x(2);sig_0(1);sig_0(2);sig_0(3)];
      
    sol = ode15s(@(t,y) ac.non_flat_model(t, y), tspan, y0, options);
    t = sol.x'; y = sol.y';
    
    % get number of complete rounds
    n_max = floor(t(end)/ac.tf);
    
    % nominal trajectory over a period
    y_nom = zeros(length(t), 6);
    for i = 1:length(t)
        sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
        ac = ac.get_xu(sig);
        y_nom(i,:) = [ac.x(1), ac.x(3), ac.x(2), sig(1), sig(2), sig(3)];
    end
   
    dev_data = y - y_nom;
end