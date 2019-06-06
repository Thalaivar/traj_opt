function [c,ceq,dc,dceq] = Cfun(Z,ac,N,D)

    VR = Z(end-1);
    tfin = Z(end);
    chi0 = Z(end-2);
    t = ((tfin/N)*(1:N))';
    
    X = zeros(N,8);
    for i = 1:8
        X(:,i) = Z((i-1)*N+1:i*N);
    end
    
    F = zeros(N,6);
    for j = 1:N
        Xval = X(j,:) + [0,chi0*t(j),zeros(1,6)];
        F(j,:) = ac.stateDervs([Xval,0],VR);
    end
    F(:,2) = F(:,2)-chi0*ones(N,1);
    
    ceq = reshape((2*pi/tfin)*D*X(:,1:6)-F,[6*N,1]);
    ceq(end+1) = chi0*tfin+2*pi;
    
    dCL = (2*pi/tfin)*D*X(:,7);
    ddCL = (2*pi/tfin)*D*dCL;
    dmu = (2*pi/tfin)*D*X(:,8);
    ddmu = (2*pi/tfin)*D*dmu;
    
    c = [];
    
    c(end+1:end+N) = X(:,6) + 0.5*ac.b*abs(sin(X(:,8))); % wing tip clearance
    
    c(end+1:end+N) = dCL - 0.4;
    c(end+1:end+N) = -dCL - 0.4;
%     c(end+1:end+N) = ddCL - 0.4;
%     c(end+1:end+N) = -ddCL - 0.4;    
    c(end+1:end+N) = dmu - (20*pi/180);
    c(end+1:end+N) = -dmu - (20*pi/180);
    c(end+1:end+N) = ddmu - (20*pi/180);
    c(end+1:end+N) = -ddmu - (20*pi/180);
    
   if nargout > 2 % gradient of the constraints
        dc = [];
        dceq = [];
    end    
end