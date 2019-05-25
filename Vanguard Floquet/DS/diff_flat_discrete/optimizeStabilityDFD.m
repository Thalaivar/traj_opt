function f = optimizeStabilityDFD(X, trajParams, windShear, objType)
    if strcmp(windShear, 'same')
        N = (length(X)-1)/3;
        VR = trajParams.VR;
    else
        N = (length(X)-2)/3;
        VR = X(3*N+2);
        trajParams.VR = VR;
    end
    T = X(3*N+1); trajParams.T = T;
    
    if strcmp(objType, 'stability')
        diffFac = 2*pi/T;
        D = trajParams.D; DD = trajParams.DD;
        x = [X(1:N), X(N+1:2*N), X(2*N+1:3*N)];
        dx = diffFac*D*[x(:,1), x(:,2), x(:,3)];
        ddx = (diffFac^2)*DD*[x(:,1), x(:,2), x(:,3)];

        p = trajParams.p;
        state = zeros(N,3); U = zeros(N,3);
        for i = 1:N
            [xx, uu] = diffFlatModel([x(i,:);dx(i,:);ddx(i,:)], p, VR);
            state(i,:) = xx; U(i,:) = uu;
        end
        
        state = [state, x];
%         FE = spectralMethod(trajParams, state, U, N);
        FE = real(freidmannMethod(trajParams, state, U));
        f = max(FE(find(abs(FE) > 1e-8)));
%         f = 0;
%         for i = 1:length(FE)
%             f = f + atan(10*FE(i));
%         end
        
    elseif strcmp(objType, 'VR')
        f = VR;
    end
end

function FE = spectralMethod(trajParams, x, u, N)
    t = trajParams.T*trajParams.t/(2*pi);
    d = 4; Dmat = zeros(d*N); Mmat = zeros(d*N); D = trajParams.D;
    
    for i = 1:N
        % get state and control vector for jacobian
        U = [u(i,:), trajParams.VR];
        A = sysModel(t(i), trajParams.p, x(i,:), U);
        for j = 1:d
           for k = 1:d
                Mmat(i+(j-1)*N,(k-1)*N+i) = A(j,k);
           end
           Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/trajParams.T);
        end
    end

    [eigVec,eigE] = eig(Dmat-Mmat);
    eigE = -1*diag(eigE);
    FE = identify_floquet(eigE, N, trajParams.T);
end

function FE = freidmannMethod(trajParams, x, u)
    d = 4;
    t = trajParams.T*trajParams.freidmannGrid/(2*pi); 
    h = t(2) - t(1);
    FTM = eye(d);
    for i = 1:length(t)-2
        K = friedman_spectral_K(h, (t(end) - i*h),  x, u, trajParams);
        FTM = FTM*K;
    end
    FM = eig(FTM);
    FE = log(FM)/trajParams.T;
end

function K = friedman_spectral_K(h, t, X, U, trajParams)
    d = 4;
    tt  = trajParams.T*trajParams.t/(2*pi);
    Cl  = interp_sinc(tt, U(:,1), t);
    mu  = interp_sinc(tt, U(:,2), t);
    V   = interp_sinc(tt, X(:,1), t);
    chi = interp_sinc(tt, X(:,2), t);
    gam = interp_sinc(tt, X(:,3), t);
    x   = interp_sinc(tt, X(:,4), t);
    y   = interp_sinc(tt, X(:,5), t);
    z   = interp_sinc(tt, X(:,6), t);
    Z = [V, chi, gam, x, y, z];
    U = [Cl, mu, 0, trajParams.VR];
    A_psi = sysModel(t, trajParams.p, Z, U);
    A_psi_h_by_2 = sysModel(t+0.5*h, trajParams.p, Z, U);
    A_psi_h = sysModel(t+h, trajParams.p, Z, U);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end

function A = sysModel(t, p, Z, U)
    prm.m = 4.5;
    prm.S = 0.473;
    prm.CD0 = 0.0173;
    prm.CD1 = -0.0337;
    prm.CD2 = 0.0517;
    prm.p_exp = p;
    % evaluate jacobian
    A = JacEval(Z, U, prm);
    A = [A(1:3,1:3), A(1:3,6); A(6,1:3), A(6,6)];
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