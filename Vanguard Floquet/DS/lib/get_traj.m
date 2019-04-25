function sigma = get_traj(t, tf, coeffs, N)
    ax = coeffs(:,1); ay = coeffs(:,2); az = coeffs(:,3);

    x = ax(1); y = ay(1); z = az(1);
    xdot = 0; ydot = 0; zdot = 0; xddot = 0; yddot = 0; zddot = 0;
    for i = 2:N+1
        x = x + ax(i)*cos(2*pi*(i-1)*t/tf) + ax(i+N)*sin(2*pi*(i-1)*t/tf);
        y = y + ay(i)*cos(2*pi*(i-1)*t/tf) + ay(i+N)*sin(2*pi*(i-1)*t/tf);
        z = z + az(i)*cos(2*pi*(i-1)*t/tf) + az(i+N)*sin(2*pi*(i-1)*t/tf);
        xdot = xdot - ax(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ax(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
        ydot = ydot - ay(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ay(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
        zdot = zdot - az(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + az(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
        xddot = xddot - ax(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ax(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
        yddot = yddot - ay(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ay(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
        zddot = zddot - az(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - az(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
    end

    sigma = [x,y,z,xdot,ydot,zdot,xddot,yddot,zddot];
end
