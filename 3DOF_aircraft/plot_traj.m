function plot_traj(a, eta, tf)
    t = linspace(0, tf, 1000);
    a_h = a(:,1); a_x = a(:,2); a_y = a(:,3);
    eta_h = eta(:,1); eta_x = eta(:,2); eta_y = eta(:,3);
    ph = [a_h,[eta_h;0]]; px = [a_x,[eta_x;0]]; py = [a_y,[eta_y;0]];
    
    h = zeros(length(t),1); x = zeros(length(t),1); y = zeros(length(t),1);
    for i = 1:length(t)
        [h(i), x(i), y(i), ~, ~, ~, ~, ~, ~] = get_traj(t(i), ph, px, py, tf);
    end
    
    
    plot3(x, y, h, 'k', 'Linewidth', 2.0); 
    xlabel("x (m)"); ylabel("y (m)"); zlabel("z (m)");
    grid on
    
end