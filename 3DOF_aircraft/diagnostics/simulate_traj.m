addpath('../trajectory')
addpath('../constraint_funcs')
addpath('../floquet')
addpath('../solutions')
addpath('..')

%load('../solutions/trajectory_opt/expo_O.mat')

% time evolution
tspan = [0,200*ac.tf];
sig_0 = get_traj(tspan(1), ac.tf, ac.coeffs, ac.N);
ac = ac.get_xu(sig_0);
y0 = [ac.x(1);ac.x(3);ac.x(2);sig_0(1);sig_0(2);sig_0(3)];
options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11, 'Events', @(t,y) z_event(t,y));
sol = ode15s(@(t,y) ac.non_flat_model(t, y), tspan, y0, options);
t1 = sol.x'; y1 = sol.y';
options = odeset('Events', @(t,y) z_event(t,y));
sol = ode15s(@(t,y) ac.non_flat_model(t, y), tspan, y0, options);
t2 = sol.x'; y2 = sol.y';

% get number of complete rounds
n_max_1 = floor(t1/ac.tf); n_max_2 = floor(t2/ac.tf); 
% nominal trajectory
n_points = 1500;
t_nom = linspace(ac.tf/(n_points-1), ac.tf, n_points);
t_ref = zeros((n_max_1+1)*n_points+1,1);
nom_traj = zeros(6, length(t_ref));
for i = 1:length(t_ref){
    
    }
end


% 
% % interpolate to uniform grid
% n_max = floor(t(end)/ac.tf);
% t_uni = linspace(0, n_max*ac.tf, n_max*1000 + n_max + 1);
% y_uni = interp1(t, y, t_uni);
% 
% % get values for first round
% [~,i] = min(abs(t_uni - ac.tf));
% t_nom = t_uni(2:i-1); y_nom = y_uni(2:i-1,:);
% 
% % deviation of consequent rounds from first round
% for i = 1:n_max
%     [~,i1] = min(abs(t_uni - (i-1)*ac.tf));
%     [~,i2] = min(abs(t_uni - i*ac.tf));    
%     y_uni(i1+1:i2-1,:) = y_uni(i1+1:i2-1,:) - y_nom;
%     
%     y_uni(i1,:) = y_uni(i1,:) - y_nom(1,:);
% end
% y_uni(i2,:) = y_uni(i2,:) - y_nom(1,:);
% 
% % deviation of last maneouvre
% [~,i_f] = min(abs(t-n_max*ac.tf));
% tf = linspace(n_max*ac.tf, t(end), 1000);
% yf = interp1(t(i_f:end), y(i_f:end,:), tf);
% tf_nom = linspace(0, t(end)- n_max*ac.tf, 1000);
% yf_nom = interp1(t_nom, y_nom, tf_nom);
% nf = length(tf);
% yf = yf - yf_nom;
% y_uni(end+1:end+length(tf),:) = yf;
% t_uni(end+1:end+length(tf)) = tf;
% 
% plot(t_uni, y_uni(:,1))
rmpath('../solutions')
rmpath('../trajectory')
rmpath('../constraint_funcs')
rmpath('../floquet')
rmpath('..')