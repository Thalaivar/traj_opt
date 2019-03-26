addpath('../../')
addpath('../../floquet')
addpath('../../trajectory')

% load stable trajectory
load('D:\ranjth_mohan\traj_opt\3DOF_aircraft\solutions\stability_opt\exponential\EOS10_10_50.mat')

% generate state time history
t = linspace(0, ac.tf, 1000);
state = zeros(length(t), 6);
pos_dot = zeros(length(t), 3);  % time histories of xdot, ydot, zdot
for i = 1:length(t)
    sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
    ac = ac.get_xu(sig);
    state(i,:) = [ac.x(1), ac.x(3), ac.x(2), sig(1), sig(2), sig(3)];
    pos_dot(i,:) = [sig(4), sig(5), sig(6)];
end

E_a_stable = zeros(1, length(t));
E_i_stable = zeros(1, length(t));
for i = 1:length(t)
    E_a_stable(i) = 0.5*ac.m*(state(i,1)^2) + ac.m*ac.g*(-state(i,6));
    E_i_stable(i) = 0.5*ac.m*(pos_dot(i,1)^2 + pos_dot(i,2)^2 + pos_dot(i,3)^2) + ac.m*ac.g*(-state(i,6));
end

% load unstable trajectory
load('D:\ranjth_mohan\traj_opt\3DOF_aircraft\solutions\stability_opt\exponential\EOU10_1_50.mat')

% generate state time history
t = linspace(0, ac.tf, 1000);
ustate = zeros(length(t), 6);
pos_dot = zeros(length(t), 3);  % time histories of xdot, ydot, zdot
for i = 1:length(t)
    sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
    ac = ac.get_xu(sig);
    ustate(i,:) = [ac.x(1), ac.x(3), ac.x(2), sig(1), sig(2), sig(3)];
    pos_dot(i,:) = [sig(4), sig(5), sig(6)];
end

E_a_unstable = zeros(1, length(t));
E_i_unstable = zeros(1, length(t));
for i = 1:length(t)
    E_a_unstable(i) = 0.5*ac.m*(ustate(i,1)^2) + ac.m*ac.g*(-ustate(i,6));
    E_i_unstable(i) = 0.5*ac.m*(pos_dot(i,1)^2 + pos_dot(i,2)^2 + pos_dot(i,3)^2) + ac.m*ac.g*(-ustate(i,6));
end
% subplot(311)
% plot(state(:,1), E_a_stable, 'r')
% hold on
% plot(state(:,1), E_a_unstable, 'b')
% xlabel('V');
% ylabel('$E_a$', 'Interpreter', 'latex');
% title('$E_a$ vs. V for loiter', 'Interpreter', 'latex');
% legend('Stable', 'Unstable');
% 
% subplot(312)
% scatter(state(:,2), E_a_stable, 5, 'r');
% hold on
% scatter(state(:,2), E_a_unstable, 5, 'b');
% xlabel('$\chi$', 'Interpreter', 'latex');
% ylabel('$E_a$', 'Interpreter', 'latex');
% title('$E_a$ vs. $\chi$ for loiter', 'Interpreter', 'latex');
% legend('Stable', 'Unstable');
% 
% subplot(313)
% plot(state(:,3), E_a_stable, 'r')
% hold on
% plot(state(:,3), E_a_unstable, 'b')
% xlabel('$\gamma$', 'Interpreter', 'latex');
% ylabel('$E_a$', 'Interpreter', 'latex');
% title('$E_a$ vs. $\gamma$ for loiter', 'Interpreter', 'latex');
% legend('Stable', 'Unstable');

% subplot(311)
% plot(state(:,1), E_i_stable, 'r')
% hold on
% plot(state(:,1), E_i_unstable, 'b')
% xlabel('V');
% ylabel('$E_i$', 'Interpreter', 'latex');
% title('$E_i$ vs. V for loiter', 'Interpreter', 'latex');
% legend('Stable', 'Unstable');
% 
% subplot(312)
% scatter(state(:,2), E_i_stable, 5, 'r');
% hold on
% scatter(state(:,2), E_i_unstable, 5, 'b');
% xlabel('$\chi$', 'Interpreter', 'latex');
% ylabel('$E_i$', 'Interpreter', 'latex');
% title('$E_i$ vs. $\chi$ for loiter', 'Interpreter', 'latex');
% legend('Stable', 'Unstable');
% 
% subplot(313)
% plot(state(:,3), E_i_stable, 'r')
% hold on
% plot(state(:,3), E_i_unstable, 'b')
% xlabel('$\gamma$', 'Interpreter', 'latex');
% ylabel('$E_i$', 'Interpreter', 'latex');
% title('$E_i$ vs. $\gamma$ for loiter', 'Interpreter', 'latex');
% legend('Stable', 'Unstable');

rmpath('../../')
rmpath('../../floquet')
rmpath('../../trajectory')