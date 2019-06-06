function [c,ceq,dc,dceq] = Cfun(Z,ac,N)

    VR = Z(end-1);
    tfin = Z(end);
    Psi0 = Z(end-2);
    t = ((tfin/N)*(1:N))';
    fac = 2*pi/tfin; fac2 = fac*fac;

    x = Z(1:N);  
    y = Z(N+1:2*N); 
    z = Z(2*N+1:3*N); z(z>0)=0;
    Phi = Z(3*N+1:4*N);  
    Thet = Z(4*N+1:5*N); 
    Psi1 = Z(5*N+1:6*N); 
    Psi = Psi1 + Psi0*t;
    
    [dx,ddx] = dftDerv(x,tfin);
    [dy,ddy] = dftDerv(y,tfin);
    [dz,ddz] = dftDerv(z,tfin);
    [dPhi,ddPhi] = dftDerv(Phi,tfin);
    [dThet,ddThet] = dftDerv(Thet,tfin);
    [dPsi1,ddPsi1] = dftDerv(Psi1,tfin);
    dPsi = dPsi1 + Psi0;
    ddPsi = ddPsi1;

    CTx(N,1) = 0; de(N,1) = 0; df(N,1) = 0; 
    VT(N,1) = 0; bet(N,1) = 0; CL(N,1) = 0;
    da(N,1) = 0; dr(N,1) = 0; CTy(N,1) = 0; 
    p(N,1) = 0; aoa(N,1) = 0;
    
    for j = 1:N
        
        Wx = VR*(-z(j)).^ac.p_exp;
        Wxz = (ac.p_exp*VR)*((-z(j)).^ac.p_exp)./z(j);
        wind = [Wx,Wxz];
        
        X = evalZ_6DOF([Phi(j),Thet(j),Psi(j),x(j),y(j),z(j)],[dPhi(j),dThet(j),dPsi(j),dx(j),dy(j),dz(j);...
        ddPhi(j),ddThet(j),ddPsi(j),ddx(j),ddy(j),ddz(j)],wind,ac);
        u = X(1);
        v = X(2);
        w = X(3);
        p(j) = X(4);
        q = X(5);
        df(j) = X(7);
        da(j) = X(8);
        de(j) = X(9);
        dr(j) = X(10);
        CTx(j) = X(11);
        CTy(j) = X(12);

        VT(j) = sqrt(u^2 + v^2 + w^2);
        aoa(j) = atan(w/u);
        bet(j) = asin(v/VT(j));

        CL(j) = ( ac.CL0 + ac.CLalf*aoa(j) + ac.CLq*(0.5*q*ac.c/VT(j)) + ac.CLdf*df(j) + ac.CLde*de(j));        
        
    end  
    
    [d_de,dd_de] = dftDerv(de,tfin);
    [d_df,dd_df] = dftDerv(df,tfin);
    [d_da,dd_da] = dftDerv(da,tfin);
    [d_dr,dd_dr] = dftDerv(dr,tfin);

%     inequality constraints
    c = [];
    c(end+1:end+N) = z + 0.5*ac.b*abs(sin(Phi)); % wing tip clearance
    c(end+1:end+N) = CL - 1.17;
    c(end+1:end+N) = -CL - (0.2);
%     c(end+1:end+N) = CTx - 1e-4;
%     c(end+1:end+N) = -CTx - 1e-4;
    c(end+1:end+N) = VT - 80;
    c(end+1:end+N) = -VT + 10;
    c(end+1:end+N) = bet - 2.0*pi/180;              % beta
    c(end+1:end+N) = -bet - 2.0*pi/180;             % beta
    c(end+1:end+N) = df - (45*pi/180);              % df
    c(end+1:end+N) = -df - (45*pi/180);             % df
    c(end+1:end+N) = da - (45*pi/180);              % da
    c(end+1:end+N) = -da - (45*pi/180);             % da    
    c(end+1:end+N) = de - (45*pi/180);              % de
    c(end+1:end+N) = -de - (45*pi/180);             % de    
    c(end+1:end+N) = dr - (45*pi/180);              % dr
    c(end+1:end+N) = -dr - (45*pi/180);             % dr 
    
    
    c(end+1:end+N) = d_df - pi/4;
    c(end+1:end+N) = -d_df - pi/4;
    c(end+1:end+N) = d_de - pi/4;
    c(end+1:end+N) = -d_de - pi/4;
    c(end+1:end+N) = d_da - pi/5;
    c(end+1:end+N) = -d_da - pi/5;
    c(end+1:end+N) = d_dr - pi/5;
    c(end+1:end+N) = -d_dr - pi/5;    

    c(end+1:end+N) = dd_df - pi/4;
    c(end+1:end+N) = -dd_df - pi/4;
    c(end+1:end+N) = dd_de - pi/4;
    c(end+1:end+N) = -dd_de - pi/4;
    c(end+1:end+N) = dd_da - pi/15;
    c(end+1:end+N) = -dd_da - pi/15;
    c(end+1:end+N) = dd_dr - pi/9;
    c(end+1:end+N) = -dd_dr - pi/9;
    
%     c(end+1:end+N) = dd_df - pi/30;
%     c(end+1:end+N) = -dd_df - pi/30;
%     c(end+1:end+N) = dd_de - pi/30;
%     c(end+1:end+N) = -dd_de - pi/30;
%     c(end+1:end+N) = dd_da - pi/30;
%     c(end+1:end+N) = -dd_da - pi/30;
%     c(end+1:end+N) = dd_dr - pi/30;
%     c(end+1:end+N) = -dd_dr - pi/30;
    
%     c(end+1:end+N) = 0.5*1.225*(VT.^2).*CL*ac.S/(ac.m*9.806) - 2; % load factor
%     c(end+1:end+N) = CTy - 1e-4;
%     c(end+1:end+N) = -CTy - 1e-4;
    c(end+1:end+N) = 0.5*p.*ac.b./VT - 0.1;
    c(end+1:end+N) = -0.5*p.*ac.b./VT - 0.1;
    c(end+1:end+N) = aoa - 15*pi/180;
    c(end+1:end+N) = -aoa - 15*pi/180;
    
%     equality constraint
    ceq = [];
    ceq(end+1:end+N) = CTx;
    ceq(end+1:end+N) = CTy;    
%     ceq(end+1:end+N) = CTx.^2;
%     ceq(end+1:end+N) = CTy.^2;
%     ceq(end+1) = 1e-3*norm(CTx);
%     ceq(end+1) = 1e-3*norm(CTy);

    ceq(end+1) = Psi0*tfin+2*pi; % for loiter
%     ceq(end+1) = Psi0; % for eight
    
    if ~isreal(ceq)||~isreal(c)
    end
    if max(isnan(c))||max(isnan(ceq))
    end
    
%     c(isnan(c)) = 1000;
%     c(isinf(c)) = -1000;
%     ceq(isnan(ceq)) = 1000;
%     ceq(isinf(ceq)) = -1000;
    
%     figure
%     subplot(2,1,1)
%     plot(ceq)
%     title('C_{eq}')    
%     subplot(2,1,2)
%     plot(c)
%     title('C');
    
    if nargout > 2 % gradient of the constraints
        dc = [];
        dceq = [];
    end    
end