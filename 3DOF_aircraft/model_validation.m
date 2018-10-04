u = [Cl, nu, T];
par = [rho, Cd0, Cd1, Cd2, S, m, g, VR, hR];

Y0 = [10; 0; 0; 20; 0; 0];
[t, y] = ode45(@(t, y) aircraft_model(t, y, u, par), cheb_t, Y0);

h_time_prop = y_time_prop(:,4);
hdot_time_prop = (-2/tf)*D*h_time_prop;
hddot_time_prop = (-2/tf)*D*hdot_time_prop;
x_time_prop = y_time_prop(:,5);
xdot_time_prop = (-2/tf)*D*x_time_prop;
xddot_time_prop = (-2/tf)*D*xdot_time_prop;
Y_time_prop = y_time_prop(:,5);
Ydot_time_prop = (-2/tf)*D*Y_time_prop;
Yddot_time_prop = (-2/tf)*D*Ydot_time_prop;
for i = 1:N+1
z = [h_time_prop(i), x_time_prop(i), Y_time_prop(i), hdot_time_prop(i), xdot_time_prop(i), Ydot_time_prop(i), hddot_time_prop(i), xddot_time_prop(i), Yddot_time_prop(i)];
end

N = n_cheb_points;
state_df_model = zeros(N+1,7);
model_par = [m, rho, S, g, Cd0, Cd1, Cd2];
for i = 1:N+1
z = [h_time_prop(i), x_time_prop(i), Y_time_prop(i), hdot_time_prop(i), xdot_time_prop(i), Ydot_time_prop(i), hddot_time_prop(i), xddot_time_prop(i), Yddot_time_prop(i)];
[V, gamma, psi, nu, Cl, T, CT] = DF_aircraft_model(z, VR, model_par);
state_df_model(i,:) = [V, gamma, psi, nu, Cl, T, CT];
end 
