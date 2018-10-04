function c = ineq_constr(z, lim, wind_par, model_par)
    % z should be of the form:
    %       [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot] = z;
    
    % limits should be of the form:
    %       limits = [Clmax, Vmax, nu_min, nu_max, Tmin, Tmax, hmin]
    
    % wind_par should be of the form:
    %       wind_par = [VR, hR]
    
    % model_par should be of the form:
    %       model_par = [m, rho, S, g, b, Cd0, Cd1, Cd2]
    
    Clmax = lim(1); Vmax = lim(2); nu_min = lim(3); nu_max = lim(4);
    CTmin = lim(6); CTmax = lim(5); hmin = lim(7);
    
    b = model_par(end); % wingspan
    h = z(1);
    [V, ~, ~, nu, Cl, ~, CT] = DF_aircraft_model(z, wind_par, model_par);
    c(1,1) = Cl - Clmax;
    c(2,1) = V - Vmax;
    c(3,1) = nu_min - nu;
    c(4,1) = nu - nu_max;
    c(5,1) = CTmin - CT;
    c(6,1) = CT - CTmax;
    c(7,1) = hmin - h + 0.5*b*sin(nu);
    c(8,1) = hmin - h - 0.5*b*sin(nu);
    
end