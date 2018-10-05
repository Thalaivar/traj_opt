function [c, ceq] = constFun(x, limits, model_par, N)
    % limits is of the form:
    %       limits = [Clmax, Vmax, nu_min, nu_max, Tmin, Tmax, hmin]    
    
    % model_par is of the form:
    %       model_par = [m, rho, S, g, b, Cd0, Cd1, Cd2]
    
    % z is of the form:
    %       z = [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot]
    
    % n1 : n_coeffs; n2 : n_phase_angles;
    n1 = N+1; n2 = N;
    a_h = x(1:n1,1); a_x = x(n1+1:2*n1,1); a_y = x(2*n1 + 1:3*n1,1);
    eta_h = x(3*n1+1:3*n1+n2,1); eta_x = x((3*n1+n2+1):(3*n1+2*n2),1); eta_y = x((3*n1+2*n2+1):3*(n1+n2),1);
    ph = [a_h,[eta_h;0]]; px = [a_x,[eta_x;0]]; py = [a_y,[eta_y;0]]; VR = x(end-1,1); tf = x(end,1);
        
    wind_par = [VR, model_par(end)];
    model_par(end)
    M = 4*N;
    c = zeros(8*M,1);
    for i = 1:M+1
        j = (i-1)*8 + 1;
        t = (i-1)*tf/M;
        [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot] = get_traj(t, ph, px, py, tf);
        z = [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot];
        c(j:j+7) = ineq_constr(z, limits, wind_par, model_par);
    end
    
    [h0, x0, y0, hdot0, xdot0, ydot0, hddot0, xddot0, yddot0] = get_traj(0, ph, px, py, tf);
    z0 = [h0, x0, y0, hdot0, xdot0, ydot0, hddot0, xddot0, yddot0];
    [~, ~, psi, ~, ~, ~] = DF_aircraft_model(z0, wind_par, model_par);
    ceq = zeros(3,1);
    ceq(1,1) = psi;
    ceq(2,1) = z0(2);
    ceq(3,1) = z0(3);
end

