function [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot] = get_traj(t, ph, px, py, tf)
    
    a_h = ph(:,1); eta_h = ph(1:end-1,2);
    a_x = px(:,1); eta_x = px(1:end-1,2);
    a_y = py(:,1); eta_y = py(1:end-1,2);
    N = length(a_h);
    
    h = a_h(1); x = a_x(1); y = a_y(1);
    hdot = 0; xdot = 0; ydot = 0; 
    hddot = 0; xddot = 0; yddot = 0; 
    
    
    for i = 2:N
        h = h + a_h(i)*sin((2*pi*i*t/tf) + eta_h(i-1));
        x = x + a_x(i)*sin((2*pi*i*t/tf) + eta_x(i-1));
        y = y + a_y(i)*sin((2*pi*i*t/tf) + eta_y(i-1));
        hdot = hdot + (a_h(i)*2*pi*i/tf)*cos((2*pi*i*t/tf) + eta_h(i-1));
        xdot = xdot + (a_x(i)*2*pi*i/tf)*cos((2*pi*i*t/tf) + eta_x(i-1));
        ydot = ydot + (a_y(i)*2*pi*i/tf)*cos((2*pi*i*t/tf) + eta_y(i-1));
        hddot = hddot - a_h(i)*((2*pi*i/tf)^2)*sin((2*pi*i*t/tf) + eta_h(i-1));
        xddot = xddot - a_x(i)*((2*pi*i/tf)^2)*sin((2*pi*i*t/tf) + eta_x(i-1));
        yddot = yddot - a_y(i)*((2*pi*i/tf)^2)*sin((2*pi*i*t/tf) + eta_y(i-1));
    end
end