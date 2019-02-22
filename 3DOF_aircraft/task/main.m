addpath('..');
addpath('../solutions');
addpath('../floquet');
addpath('../trajectory');
 

D = analyze_FTM(ac);

%options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11, 'Events', @(t,y) z_event(t,y));
%data1 = simulate_traj(ac, 200, options);
% options = odeset('Events', @(t,y) z_event(t,y));
% data2 = simulate_traj(ac, 60, options);

rmpath('..');
rmpath('../solutions');
rmpath('../floquet');
rmpath('../trajectory');
