function [FTM,FE] = timeMarchMethod(trajData)
    d = 4;
    Phi=zeros(d);
    Phi0=eye(d);
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    for i=1:d
        [~,X] = ode15s(@(t,X) timeMarchModel(t, X, trajData),[0,trajData.T],Phi0(:,i),options);
        Phi(:,i)=X(end,:)';
    end
    FTM = Phi;
    FM = eig(FTM);
    FE = log(FM)/trajData.T;
end

function dX = timeMarchModel(t, X, trajData)
    A = sysModel(t, trajData);
    dX = A*X;
end

function A = sysModel(t, trajData)
    prm.m = 4.5;
    prm.S = 0.473;
    prm.CD0 = 0.0173;
    prm.CD1 = -0.0337;
    prm.CD2 = 0.0517;
    prm.p_exp = trajData.p;
    % get state and control vector for jacobian
    tt  = trajData.T*trajData.fourierGrid/(2*pi);
    V     = interp_sinc(tt, trajData.X(:,1), t);
    chi   = interp_sinc(tt, trajData.X(:,2), t);
    gamma = interp_sinc(tt, trajData.X(:,3), t);
    x     = interp_sinc(tt, trajData.X(:,4), t);
    y     = interp_sinc(tt, trajData.X(:,5), t);
    z     = interp_sinc(tt, trajData.X(:,6), t);
    CL = interp_sinc(tt, trajData.U(:,1), t);
    mu = interp_sinc(tt, trajData.U(:,2), t);
    if strcmp(trajData.type, 'diff-flat')
        CT = interp_sinc(tt, trajData.U(:,3), t);
    else
        CT = 0;
    end
    Z = [V, chi, gamma, x, y, z];
    U = [CL, mu, CT, trajData.VR];
    % evaluate jacobian
    A = JacEval(Z, U, prm);
    A = [A(1:3,1:3), A(1:3,6); A(6,1:3), A(6,6)];
end