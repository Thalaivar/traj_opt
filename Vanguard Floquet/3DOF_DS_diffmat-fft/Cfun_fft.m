function [c,ceq,dc,dceq] = Cfun_fft(X,ac,N)

    VR = X(end);
    tf = X(end-1);
%     t = linspace(tf/N,tf,N);
    fac = 2*pi/tf;
    
    x = X((1-1)*N+1:1*N);  
    y = X((2-1)*N+1:2*N); 
    z = X((3-1)*N+1:3*N); 
    
    xhat = fft(x); dxhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*xhat; % 1 is odd
    dx = fac*real(ifft(dxhat));
    ddxhat = ((1i*[0:N/2 -N/2+1:-1]').^2).*xhat; % 2 is even
    ddx = fac*fac*real(ifft(ddxhat));
    
    yhat = fft(y); dyhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*yhat; % 1 is odd
    dy = fac*real(ifft(dyhat));
    ddyhat = ((1i*[0:N/2 -N/2+1:-1]').^2).*yhat; % 2 is even
    ddy = fac*fac*real(ifft(ddyhat));    

    zhat = fft(z); dzhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*zhat; % 1 is odd
    dz = fac*real(ifft(dzhat));
    ddzhat = ((1i*[0:N/2 -N/2+1:-1]').^2).*zhat; % 2 is even
    ddz = fac*fac*real(ifft(ddzhat));
    
    V(N,1) = 0; chi(N,1) = 0; gam(N,1) = 0;
    CL(N,1) = 0; mu(N,1) = 0; CT(N,1) = 0;
    for j = 1:N
        %%% Wind model
        % exponential
        p_exp = 1;
        Wx = VR*(-z(j)).^p_exp;
        Wxz = (p_exp*VR)*((-z(j)).^p_exp)./z(j);
        wind = [Wx,Wxz];
        
        Z = evalZ_3DOF([x(j),y(j),z(j)],[dx(j),dy(j),dz(j);ddx(j),ddy(j),ddz(j)],wind,ac);
        V(j) = Z(1);
        chi(j) = Z(2);
        gam(j) = Z(3);
        CL(j) = Z(4);
        mu(j) = Z(5);
        CT(j) = Z(6);        
    end
    
    CLhat = fft(CL); dCLhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*CLhat; % 1 is odd
    dCL = fac*real(ifft(dCLhat));
    muhat = fft(mu); dmuhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*muhat; % 1 is odd
    dmu = fac*real(ifft(dmuhat));
    ddmuhat = ((1i*[0:N/2 -N/2+1:-1]').^2).*muhat; % 2 is even
    ddmu = fac*fac*real(ifft(ddmuhat));    

    % inequality constraints    
    c(11*N,1) = -1;
    c((1-1)*N+1:1*N) = z + 0.5*ac.b*abs(sin(mu)); % wing tip clearance
    c((2-1)*N+1:2*N) = CL - 1.17;
    c((3-1)*N+1:3*N) = -CL - 0.2;
%     c((4-1)*N+1:4*N) = (1e4)*CT - 1;
%     c((5-1)*N+1:5*N) = -(1e4)*CT - 1;
    c((6-1)*N+1:6*N) = V - 80;
    c((7-1)*N+1:7*N) = -V + 10;
    c((8-1)*N+1:8*N) = mu - pi/3;
    c((9-1)*N+1:9*N) = -mu - pi/3;
    c((10-1)*N+1:10*N) = gam - pi/4;
    c((11-1)*N+1:11*N) = -gam - pi/4;
    %%%% control rate constraints
    c((12-1)*N+1:12*N) = dCL - 0.5;
    c((13-1)*N+1:13*N) = -dCL - 0.5;
    c((14-1)*N+1:14*N) = dmu - (5*pi/180);
    c((15-1)*N+1:15*N) = -dmu - (5*pi/180);
    c((16-1)*N+1:16*N) = ddmu - (20*pi/180);
    c((17-1)*N+1:17*N) = -ddmu - (20*pi/180);    
    
    
%     % equality constraints
%     ceq(1) = x(1);
%     ceq(2) = y(1);
    ceq = CT;
%     ceq = [];
    
    if nargout > 2 % gradient of the constraints
      dc = [];
      dceq = [];
    end        



end