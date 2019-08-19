function [FE, FTM] = timeMarchMethodGen(trajData, jacEval, d)
    Phi=zeros(d);
    Phi0=eye(d);
    options = odeset('RelTol',1e-14,'AbsTol',1e-14);
    for i=1:d
        [~,X] = ode15s(@(t,X) timeMarchModel(t, X, trajData, jacEval),[0,trajData.T],Phi0(:,i),options);
        Phi(:,i)=X(end,:)';
    end
    FTM = Phi;
    FM = eig(FTM);
    FE = log(FM)/trajData.T;
end

function dX = timeMarchModel(t, X, trajData, jacEval)
    A = jacEval(t, trajData);
    dX = A*X;
end

