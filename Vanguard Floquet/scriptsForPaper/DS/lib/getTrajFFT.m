function [state, control] = getTrajFFT(xyz, tf, VR, p)   
    fac = 2*pi/tf;
    N = length(xyz);
    
    x = xyz(:,1); y = xyz(:,2); z = xyz(:,3); 
    xhat = fft(x); dxhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*xhat; % 1 is odd
    dx = fac*real(ifft(dxhat));
    ddxhat = ((1i*[0:N/2 -N/2+1:-1]').^2).*xhat; % 2 is even
    ddx = fac*fac*real(ifft(ddxhat));
    
    yhat = fft(y); dyhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*yhat; % 1 is odd
    dy = fac*real(ifft(dyhat));
    ddyhat = ((1i*[0:N/2 -N/2+1:-1]').^2).*yhat; % 2 is even
    ddy = fac*fac*real(ifft(ddyhat));    

    zhat = fft(z); dzhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*zhat; % 1 is odd
    dz = fac*real(ifft(dzhat));
    ddzhat = ((1i*[0:N/2 -N/2+1:-1]').^2).*zhat; % 2 is even
    ddz = fac*fac*real(ifft(ddzhat));
    
    state = zeros(N,6); control = zeros(N,3);
    for i = 1:N
        Wx = VR*(-z(i))^p;
        Wxz = (p*VR)*((-z(i))^p)/z(i);
        Y  = [x(i), y(i), z(i)];
        dY = [dx(i), dy(i), dz(i); ddx(i), ddy(i), ddz(i)];
        Z = get_xu(Y, dY, [Wx, Wxz]);
        state(i,:) = [Z(1), Z(2), Z(3), x(i), y(i), z(i)];
        control(i,:) = [Z(4), Z(5), Z(6)];
    end
end

function Z = get_xu(Y, dY, wind)
    m = 4.5;
    rho = 1.225;
    S = 0.473;
    g = 9.806;
    Cd0 = 0.0173;
    Cd1 = -0.0337;
    Cd2 = 0.0517;
    
    z = Y(3);
    zdot = dY(1,3); xdot = dY(1,1); ydot = dY(1,2); 
    zddot = dY(2,3); xddot = dY(2,1); yddot = dY(2,2);

    % wind model
    Wx = wind(1); Wxz = wind(2);

    % non flat outputs
    V = ((xdot - Wx)^2 + ydot^2 + zdot^2)^0.5;
    Vdot = (xdot*xddot - xdot*zdot*Wxz - xddot*Wx + Wx*Wxz*zdot + ydot*yddot + zdot*zddot)/V;
    gamma = -asin(zdot/V);
    gammadot = (zdot*Vdot - V*zddot)/(V*(V^2 - zdot^2)^0.5);
    chi = atan2(ydot,(xdot - Wx));
    chidot = (xdot*yddot - yddot*Wx - ydot*xddot + ydot*zdot*Wxz)/(ydot^2 + xdot^2 + Wx^2 - 2*xdot*Wx);
    nu = atan((V*cos(gamma)*chidot - Wxz*zdot*sin(chi))/(V*gammadot + g*cos(gamma) - Wxz*cos(chi)*sin(gamma)*zdot));
    Cl = (m*V*cos(gamma)*chidot - m*Wxz*zdot*sin(chi))/(0.5*rho*S*sin(nu)*V^2);
    
    % aerodynamic forces
    Cd = Cd0+ Cd1*Cl + Cd2*Cl^2;
    D = 0.5*rho*S*V^2*Cd;
    T = m*Vdot + D + m*g*sin(gamma) + m*Wxz*zdot*cos(gamma)*cos(chi);
    CT = T/(0.5*rho*S*V^2);
    Z = [V, chi, gamma, Cl, nu, CT];
end