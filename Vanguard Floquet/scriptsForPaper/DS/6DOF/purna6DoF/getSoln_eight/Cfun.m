function [c,ceq,dc,dceq] = Cfun(Z,ac,N,D)
   VR = Z(end-1);
   tfin = Z(end);
%    t = ((tfin/N)*(1:N))';
   
    % X = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z,df,da,de,dr,CTx,CTy]
    X(N,18) = 0;
    for j = 1:18
        X(:,j) = Z((j-1)*N+1:j*N);
    end
    
    F(N,12) = 0;
    for i = 1:N
        F(i,:) = ac.stateDervs(X(i,:),VR);
    end
    
    dX = (2*pi/tfin)*D*X(:,1:12);
    ceq = reshape(dX-F,[12*N,1]);
     
    df = X(:,13); [d_df,dd_df] = dftDerv(df,tfin);
    da = X(:,14); [d_da,dd_da] = dftDerv(da,tfin);
    de = X(:,15); [d_de,dd_de] = dftDerv(de,tfin);
    dr = X(:,16); [d_dr,dd_dr] = dftDerv(dr,tfin);
    
    VT = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2);
    aoa = atan(X(:,3)./X(:,1));
    bet = asin(X(:,2)./VT);
    CL = ( ac.CL0 + ac.CLalf*aoa + ac.CLq*(0.5*X(:,5)*ac.c./VT) + ac.CLdf*df + ac.CLde*de );
    
    c = [];
    c(end+1:end+N) = X(:,12) + 0.5*ac.b*abs(sin(X(:,7)));
    c(end+1:end+N) = CL - 1.17;
    c(end+1:end+N) = -CL - (0.2);
    c(end+1:end+N) = VT - 80;
    c(end+1:end+N) = -VT + 10;
    c(end+1:end+N) = bet - 2.0*pi/180;              % beta
    c(end+1:end+N) = -bet - 2.0*pi/180;             % beta

    c(end+1:end+N) = d_df - pi/4;
    c(end+1:end+N) = -d_df - pi/4;
    c(end+1:end+N) = d_de - pi/4;
    c(end+1:end+N) = -d_de - pi/4;
    c(end+1:end+N) = d_da - pi/5;
    c(end+1:end+N) = -d_da - pi/5;
    c(end+1:end+N) = d_dr - pi/5;
    c(end+1:end+N) = -d_dr - pi/5;
    
%     c(end+1:end+N) = dd_df - pi/4;
%     c(end+1:end+N) = -dd_df - pi/4;
%     c(end+1:end+N) = dd_de - pi/4;
%     c(end+1:end+N) = -dd_de - pi/4;
%     c(end+1:end+N) = dd_da - pi/5;  
%     c(end+1:end+N) = -dd_da - pi/5;
%     c(end+1:end+N) = dd_dr - pi/5;
%     c(end+1:end+N) = -dd_dr - pi/5;    

%     c(end+1:end+N) = dd_df - pi/4;
%     c(end+1:end+N) = -dd_df - pi/4;
%     c(end+1:end+N) = dd_de - pi/4;
%     c(end+1:end+N) = -dd_de - pi/4;
%     c(end+1:end+N) = dd_da - pi/15;
%     c(end+1:end+N) = -dd_da - pi/15;
%     c(end+1:end+N) = dd_dr - pi/9;
%     c(end+1:end+N) = -dd_dr - pi/9;
   

    limSet = pi/150;
    c(end+1:end+N) = dd_df - limSet;
    c(end+1:end+N) = -dd_df - limSet;
    c(end+1:end+N) = dd_de - limSet;
    c(end+1:end+N) = -dd_de - limSet;
    c(end+1:end+N) = dd_da - limSet;
    c(end+1:end+N) = -dd_da - limSet;
    c(end+1:end+N) = dd_dr - limSet;
    c(end+1:end+N) = -dd_dr - limSet;
    
%     c(end+1:end+N) = 0.5*1.225*(VT.^2).*CL*ac.S/(ac.m*9.806) - 2; % load factor
    c(end+1:end+N) = 0.5*X(:,4).*ac.b./VT - 0.1;
    c(end+1:end+N) = -0.5*X(:,4).*ac.b./VT - 0.1;
    c(end+1:end+N) = aoa - 15*pi/180;
    c(end+1:end+N) = -aoa - 15*pi/180;
    
   if nargout > 2 % gradient of the constraints
        dc = [];
        dceq = [];
    end        
end