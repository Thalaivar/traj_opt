function [c,ceq,dc,dceq] = Cfun(Z,ac,N,D)
   VR = Z(end-1);
   tfin = Z(end);
   Psi0 = Z(end-2);
   t = ((tfin/N)*(1:N))';
   
    % X = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z,df,da,de,dr,CTx,CTy]
    X(N,16) = 0;
    for j = 1:16
        X(:,j) = Z((j-1)*N+1:j*N);
    end
    X(:,9) = X(:,9) + Psi0*t;
    
    F(N,12) = 0;
    for i = 1:N
        F(i,:) = ac.stateDervs(X(i,:),VR);
    end
    
    dX = (2*pi/tfin)*D*X(:,1:12); dX(:,9) = (2*pi/tfin)*D*Z(8*N+1:9*N) + Psi0;
    ceq = reshape(dX-F,[12*N,1]);
    
    ceq(end+1) = Psi0*tfin+2*pi; % loiter
    
%     df = X(:,13); [d_df,dd_df] = dftDerv(df,tfin);
%     da = X(:,14); [d_da,dd_da] = dftDerv(da,tfin);
%     de = X(:,15); [d_de,dd_de] = dftDerv(de,tfin);
%     dr = X(:,16); [d_dr,dd_dr] = dftDerv(dr,tfin);
    
    df = X(:,13); de = X(:,15);
    
    VT = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2);
    aoa = atan(X(:,3)./X(:,1));
    bet = asin(X(:,2)./VT);
    CL = ( ac.CL0 + ac.CLalf*aoa + ac.CLq*(0.5*X(:,5)*ac.c./VT) + ac.CLdf*df + ac.CLde*de );
    
    c = [];
    c(end+1:end+N,1) = X(:,12) + 0.5*ac.b*abs(sin(X(:,7)));
    c(end+1:end+N) = CL - 1.17;
    c(end+1:end+N) = -CL - (0.2);
    c(end+1:end+N) = VT - 80;
    c(end+1:end+N) = -VT + 10;
    c(end+1:end+N) = bet - 10*pi/180;              % beta
    c(end+1:end+N) = -bet - 10*pi/180;             % beta

%     c(end+1:end+N) = d_df - pi/4;
%     c(end+1:end+N) = -d_df - pi/4;
%     c(end+1:end+N) = d_de - pi/4;
%     c(end+1:end+N) = -d_de - pi/4;
%     c(end+1:end+N) = d_da - pi/5;
%     c(end+1:end+N) = -d_da - pi/5;
%     c(end+1:end+N) = d_dr - pi/5;
%     c(end+1:end+N) = -d_dr - pi/5;
% 
%     limSet = pi/60;
%     c(end+1:end+N) = dd_df - limSet;
%     c(end+1:end+N) = -dd_df - limSet;
%     c(end+1:end+N) = dd_de - limSet;
%     c(end+1:end+N) = -dd_de - limSet;
%     c(end+1:end+N) = dd_da - limSet;
%     c(end+1:end+N) = -dd_da - limSet;
%     c(end+1:end+N) = dd_dr - limSet;
%     c(end+1:end+N) = -dd_dr - limSet;
    
%     c(end+1:end+N) = 0.5*1.225*(VT.^2).*CL*ac.S/(ac.m*9.806) - 2; % load factor
    c(end+1:end+N) = 0.5*X(:,4).*ac.b./VT - 1;
    c(end+1:end+N) = -0.5*X(:,4).*ac.b./VT - 1;
    c(end+1:end+N) = aoa - 15*pi/180;
    c(end+1:end+N) = -aoa - 15*pi/180;
    
    %      X = zeros(N,12); U = zeros(N,4);
    %     for i = 1:12
    %         j = (i-1)*N;
    %         X(:,i) = Z(j+1:j+N);
    %     end
    %     for i = 13:16
    %         j = (i-1)*N;
    %         U(:,i-12) = Z(j+1:j+N);
    %     end
    %     trajData.N = N; trajData.X = X;
    %     trajData.U = U; trajData.D = D;
    %     trajData.ac = ac;
    %         
    %     trajData.T = Z(16*N+3); trajData.VR = Z(16*N+2);
    %     [FE, ~, AM, groupSizes] = spectralMethod(trajData);
    %     
    %     if imag(FE(1))~=0
    %         c(end+1) = 1e-2 - real(FE(1) - FE(3))/abs(real(FE(1)));
    %     else
    %         c(end+1) = 1e-2 - real(FE(1) - FE(2))/abs(real(FE(1)));
    %     end
    %     c(end+1) = groupSizes(1) - (N+6);
    %     ceq(end+1) = AM(1) - 1;

   if nargout > 2 % gradient of the constraints
        dc = [];
        dceq = [];
    end        
end