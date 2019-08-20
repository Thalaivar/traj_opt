function A = jac3DoFTimeInterp(t, trajData)
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
    if(isfield(trajData, 'chiLinearTerm'))
        chi = chi + trajData.chiLinearTerm*t;
    end
    gamma = interp_sinc(tt, trajData.X(:,3), t);
    x     = interp_sinc(tt, trajData.X(:,4), t);
    y     = interp_sinc(tt, trajData.X(:,5), t);
    z     = interp_sinc(tt, trajData.X(:,6), t);
    CL = interp_sinc(tt, trajData.U(:,1), t);
    mu = interp_sinc(tt, trajData.U(:,2), t);
    CT = 0;
%     if strcmp(trajData.type, 'diff-flat')
%         CT = interp_sinc(tt, trajData.U(:,3), t);
%     else
%         CT = 0;
%     end
    Z = [V, chi, gamma, x, y, z];
    U = [CL, mu, CT, trajData.VR];
    % evaluate jacobian
    A = JacEval(Z, U, prm);
    if(trajData.p == 0.25)
        A = [A(1:3,1:3), A(1:3,6); A(6,1:3), A(6,6)];
    elseif(trajData.p == 1)
        A = A(1:3, 1:3);
    end
end