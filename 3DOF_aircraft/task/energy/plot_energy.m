function plot_energy(max_FM, E, state, type)    
    if strcmp(type, 'air-relative')
        E_label = '$E_a$';
    else
        E_label = '$E_i$';
    end
    
    n_graphs = length(state);
    
    subplot(311)
    hold on
    for i = 1:n_graphs
        plot(state{i}(:,1), E(i,:), 'DisplayName', string(max_FM(i)))
    end
    xlabel('V');
    ylabel(E_label, 'Interpreter', 'latex');
    title(strcat(E_label, ' vs. V for loiter'), 'Interpreter', 'latex');
    legend show;
    
    subplot(312)
    hold on
    for i = 1:n_graphs
        scatter(state{i}(:,2), E(i,:), 5, 'DisplayName', string(max_FM(i)));
    end
    xlabel('$\chi$', 'Interpreter', 'latex');
    ylabel(E_label, 'Interpreter', 'latex');
    title(strcat(E_label, ' vs. $\chi$ for loiter'), 'Interpreter', 'latex');
    legend show;
    
    subplot(313)
    hold on
    for i = 1:n_graphs    
        plot(state{i}(:,3), E(i,:), 'DisplayName', string(max_FM(i)));
    end
    xlabel('$\gamma$', 'Interpreter', 'latex');
    ylabel(E_label, 'Interpreter', 'latex');
    title(strcat(E_label, ' vs. $\gamma$ for loiter'), 'Interpreter', 'latex');
    legend show;
end