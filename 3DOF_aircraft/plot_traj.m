function plot_traj(coeffs, tf, N)
    t = linspace(0, tf, 1000);
    z = zeros(1,length(t)); x = zeros(1,length(t)); y = zeros(1,length(t));
    for i = 1:length(t)
        sigma = aircraft.get_traj(t(i), tf, coeffs, N);
        x(i) = sigma(1); y(i) = sigma(2); z(i) = sigma(3);
    end
    plot3(x, -y, -z, 'k', 'Linewidth', 2.0); 
    xlabel("x (m)"); ylabel("y (m)"); zlabel("z (m)");
    grid on
    
end