function plot_lam_vs_E(E, lam, max_FM)
    n_traj = size(E); n_traj = n_traj(1);
    leg = cell(1, n_traj);
    for i = 1:n_traj
        mean_E = mean(E(i,:));
        scatter(mean_E*ones(4,1), real(lam(:,i)), 20, 'filled')
        hold on
        leg{i} = num2str(max_FM(i));
    end
    xlabel('E');
    ylabel('$real(\frac{log(\lambda)}{t_f})$', 'Interpreter', 'latex');
    title('Air-relative energy vs. Real part of FE');
    legend(leg);
end