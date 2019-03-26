addpath('../../');
addpath('../../solutions');
addpath('../../floquet');
addpath('../../trajectory');
 
% for i=1:4
%     div_tol = [1/(10^i), 1e-1];
%     D = analyze_FTM(ac, div_tol)
% end

D = analyze_FTM(ac, [1e-2, 1e-2], 'expo')

% options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11, 'Events', @(t,y) z_event(t,y));
% data1 = simulate_traj(ac, 1000, options);
% options = odeset('Events', @(t,y) z_event(t,y));
% data2 = simulate_traj(ac, 60, options);

rmpath('../../');
rmpath('../../solutions');
rmpath('../../floquet');
rmpath('../../trajectory');
