function [c, ceq] = constFun(x, limits, model_par, N)
    % limits is of the form:
    %       limits = [Clmax, Vmax, nu_min, nu_max, Tmin, Tmax, hmin]    
    
    % model_par is of the form:
    %       model_par = [m, rho, S, g, b, Cd0, Cd1, Cd2]
    
    % z is of the form:
    %       z = [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot]
    
    a_h = x(1:N+1,1); a_x = x(N+2:2*(N+1),1); a_y = x(2*N+3:3*(N+1),1);
    eta_h = x(3*N+4:4*(N+1),1); eta_x = x(4*N+5:5*(N+1),1); eta_y = x(5*N+6:6*(N+1),1);
    ph = [a_h,eta_h]; px = [a_x,eta_x]; py = [a_y,eta_y]; VR = x(end-1,1); tf = x(end,1);
    
    wind_par = VR;
    
    M = 5*N;
    c = zeros(8*M,1);
    for i = 1:8:M+1
        t = (i-1)*tf/M;
        [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot] = get_traj(t, ph, px, py, tf);
        z = [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot];
        c(i:i+7) = ineq_constr(z, limits, wind_par, model_par); 
    end
    
    [h0, x0, y0, hdot0, xdot0, ydot0, hddot0, xddot0, yddot0] = get_traj(0, ph, px, py, tf);
    z0 = [h0, x0, y0, hdot0, xdot0, ydot0, hddot0, xddot0, yddot0];
    [~, ~, psi, ~, ~, ~] = DF_aircraft_model(z0, wind_par, model_par);
    ceq = zeros(3,1);
    ceq(1,1) = psi;
    ceq(2,1) = z0(2);
    ceq(3,1) = z0(3);
end

