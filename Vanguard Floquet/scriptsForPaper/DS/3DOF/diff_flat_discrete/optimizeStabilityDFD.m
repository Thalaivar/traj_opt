function f = optimizeStabilityDFD(X, trajData, windShear, objType)
    X = X';
    if strcmp(windShear, 'same')
        N = (length(X)-1)/3;
        VR = trajData.VR;
    else
        N = (length(X)-2)/3;
        VR = X(3*N+2);
        trajData.VR = VR;
    end
    T = X(3*N+1); trajData.T = T;
    trajData.N = N;
    
    if strcmp(objType, 'stability')
        diffFac = 2*pi/T;
        D = trajData.D; DD = trajData.DD;
        x = [X(1:N), X(N+1:2*N), X(2*N+1:3*N)];
        dx = diffFac*D*[x(:,1), x(:,2), x(:,3)];
        ddx = (diffFac^2)*DD*[x(:,1), x(:,2), x(:,3)];

        p = trajData.p;
        state = zeros(N,3); U = zeros(N,3);
        for i = 1:N
            [xx, uu] = diffFlatModel([x(i,:);dx(i,:);ddx(i,:)], p, VR);
            state(i,:) = xx; U(i,:) = uu;
        end
        
        x = [state, x];
        trajData.X = x; trajData.U = U;
        [~,f] = spectralMethod(trajData, true);
%         FE = real(freidmannMethod(trajParams, state, U));
%         f = max(FE(find(abs(FE) > 1e-8)));
%         f = 0;
%         for i = 1:length(FE)
%             f = f + atan(10*FE(i));
%         end
        
    elseif strcmp(objType, 'VR')
        f = VR;
    end
end

function [x, u] = diffFlatModel(X, p, VR)
    g = 9.81;
    m = 4.5;
    rho = 1.225;
    S = 0.473;
    g = 9.806;
    Cd0 = 0.0173;
    Cd1 = -0.0337;
    Cd2 = 0.0517;
    
    z = X(1,3);
    zdot = X(2,3); xdot = X(2,1); ydot = X(2,2); 
    zddot = X(3,3); xddot = X(3,1); yddot = X(3,2);

    % wind model
    Wx = VR*(-z)^p;
    Wxz = (p*VR)*((-z)^p)/z;

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
    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    D = 0.5*rho*S*V^2*Cd;
    T = m*Vdot + D + m*g*sin(gamma) + m*Wxz*zdot*cos(gamma)*cos(chi);
    CT = T/(0.5*rho*S*V^2);
    
    x = [V, chi, gamma]; u = [Cl, nu, CT];
end