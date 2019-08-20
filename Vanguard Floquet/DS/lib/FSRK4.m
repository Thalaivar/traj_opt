function [FE, FTM] = FSRK4(trajData, FSGrid, jacEval, d)
    t = trajData.T*FSGrid/(2*pi);
    h = t(2) - t(1);
    FTM = eye(d);
    for i = 1:length(t)-2
        K = spectralK(h, (t(end) - i*h), d, trajData, jacEval);
        FTM = FTM*K;
    end
    FM = eig(FTM);
    FE = log(FM)/trajData.T;
end

function K = spectralK(h, t, d, trajData, jacEval)
    A_psi = jacEval(t, trajData);
    A_psi_h_by_2 = jacEval(t+0.5*h,trajData);
    A_psi_h = jacEval(t+h, trajData);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end