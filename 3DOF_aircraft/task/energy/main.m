 addpath('../../')
addpath('../../floquet')
addpath('../../trajectory')

clearvars

traj = {'../../unstable_close_to_1.mat', '../../stable_close_to_1.mat', '../../EOS1_01_50.mat', '../../EOS01_01_50.mat', '../../EOS5_01_50.mat', '../../EOS10_01_50.mat', '../../EOS15_01_50.mat', '../../EOU5_01_50.mat', '../../EOU10_01_50.mat', '../../EOU15_01_50.mat'};
n_traj = length(traj); n_points = 1000;
E_a = zeros(n_traj, n_points);  
E_i = zeros(n_traj, n_points);
state = cell(1, n_traj);
FE = zeros(4, n_traj);
max_FM = zeros(1, n_traj);

for i = 1:n_traj
    load(traj{i});
    [E_a_temp, E_i_temp, state_temp] = analyze_E(ac, n_points);
    E_a(i,:) = E_a_temp;
    E_i(i,:) = E_i_temp;
    state{i} = state_temp;
    FTM = get_FTM(ac, 'friedmann'); FTM = [FTM(1:3,1:3),FTM(1:3,6);FTM(6,1:3),FTM(6,6)];
    D1 = eig(FTM);
    D = log(D1)/ac.tf;
    FE(:,i) = D;
    max_FM(i) = max(abs(D1));
    
    clear('ac','p','sol');
end

plot_lam_vs_E(E_a, FE, max_FM);
    
rmpath('../../')
rmpath('../../floquet')
rmpath('../../trajectory')