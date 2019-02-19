addpath('..');
addpath('../solutions');
addpath('../floquet');
addpath('../trajectory');
 
options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11, 'Events', @(t,y) z_event(t,y));
[t1, dev_data1] = simulate_traj(ac, 200, options);
options = odeset('Events', @(t,y) z_event(t,y));
[t2, dev_data2] = simulate_traj(ac, 60, options);

rmpath('..');
rmpath('../solutions');
rmpath('../floquet');
rmpath('../trajectory');
