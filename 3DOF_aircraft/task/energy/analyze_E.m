function [E_a, E_i, state] = analyze_E(ac, n_points)
    % generate state time history
    t = linspace(0, ac.tf, n_points);
    state = zeros(length(t), 6);
    pos_dot = zeros(length(t), 3);  % time histories of xdot, ydot, zdot
    for i = 1:length(t)
        sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
        ac = ac.get_xu(sig);
        state(i,:) = [ac.x(1), ac.x(3), ac.x(2), sig(1), sig(2), sig(3)];
        pos_dot(i,:) = [sig(4), sig(5), sig(6)];
    end

    E_a = zeros(1, length(t));
    E_i = zeros(1, length(t));
    for i = 1:length(t)
        E_a(i) = 0.5*ac.m*(state(i,1)^2) + ac.m*ac.g*(-state(i,6));
        E_i(i) = 0.5*ac.m*(pos_dot(i,1)^2 + pos_dot(i,2)^2 + pos_dot(i,3)^2) + ac.m*ac.g*(-state(i,6));
    end
end