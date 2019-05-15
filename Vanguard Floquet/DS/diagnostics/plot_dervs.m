function plot_dervs(xdot, xdot_cap, N, p, T)
    if p == 1
        profile = 'linear';
    else
        profile = 'exponential';
    end
    [~,tt] = fourierdiff(N);
    tt = T*tt/(2*pi);
    % plot derivatives
    subplot(321)
    plot(tt, xdot(1:N), 'r')
    hold on
    plot(tt, xdot_cap(1:N), '--ob')
    ylabel('$V$', 'Interpreter', 'latex'); xlabel('$t$ (s)', 'Interpreter', 'latex');
    title(['$\dot{V}$ vs. $t$ for N = ', num2str(N), ' for ', profile, ' wind'], 'Interpreter', 'latex');
    legend('Actual', 'Spectral')
    subplot(322)
    plot(tt, xdot(N+1:2*N), 'r', 'LineWidth', 1.5)
    hold on
    plot(tt, xdot_cap(N+1:2*N), '--ob')
    ylabel('$\dot{\chi}$', 'Interpreter', 'latex')
    xlabel('t (s)', 'Interpreter', 'latex')
    title(['$\dot{\chi}$ vs. $t$ for N = ', num2str(N), ' for ', profile, ' wind'], 'Interpreter', 'latex');
    legend('Actual', 'Spectral')
    subplot(323)
    plot(tt, xdot(2*N+1:3*N), 'r')
    hold on
    plot(tt, xdot_cap(2*N+1:3*N), '--ob')
    legend('Actual', 'Spectral')
    ylabel('$\gamma$', 'Interpreter', 'latex')
    xlabel('t (s)', 'Interpreter', 'latex')
    title(['$\dot{\gamma}$ vs. $t$ for N = ', num2str(N), ' for ', profile, ' wind'], 'Interpreter', 'latex');
    subplot(324)
    plot(tt, xdot(3*N+1:4*N), 'r')
    hold on
    plot(tt, xdot_cap(3*N+1:4*N), '--ob')
    ylabel('$x$', 'Interpreter', 'latex')
    xlabel('t (s)', 'Interpreter', 'latex')
    legend('Actual', 'Spectral')
    title(['$\dot{x}$ vs. $t$ for N = ', num2str(N), ' for ', profile, ' wind'], 'Interpreter', 'latex');
    subplot(325)
    plot(tt, xdot(4*N+1:5*N), 'r')
    hold on
    plot(tt, xdot_cap(4*N+1:5*N), '--ob')
    legend('Actual', 'Spectral')
    ylabel('$y$', 'Interpreter', 'latex')
    xlabel('t (s)', 'Interpreter', 'latex')
    title(['$\dot{y}$ vs. $t$ for N = ', num2str(N), ' for ', profile, ' wind'], 'Interpreter', 'latex');
    subplot(326)
    plot(tt, xdot(5*N+1:6*N), 'r')
    hold on
    plot(tt, xdot_cap(5*N+1:6*N), '--ob')
    legend('Actual', 'Spectral')
    ylabel('$z$', 'Interpreter', 'latex')
    xlabel('t (s)', 'Interpreter', 'latex')
    title(['$\dot{z}$ vs. $t$ for N = ', num2str(N), ' for ', profile, ' wind'], 'Interpreter', 'latex');
end