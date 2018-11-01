function [c,ceq,dc,dceq] = Cfun(X,ac,M,N)
    
    VR = X(end);
    tf = X(end-1);
    t = linspace(0,tf,N);
    
    fourCoeff(3,2*M+1) = 0;
    len = 2*M+1;
    for i = 1:3
        fourCoeff(i,:) = X((i-1)*len+1:i*len);
    end
    
    x(N,1) = 0; y(N,1) = 0; z(N,1) = 0;
    dx(N,1) = 0; dy(N,1) = 0; dz(N,1) = 0;
    ddx(N,1) = 0; ddy(N,1) = 0; ddz(N,1) = 0;
    
    V(N,1) = 0; chi(N,1) = 0; gam(N,1) = 0;
    CL(N,1) = 0; mu(N,1) = 0; CT(N,1) = 0;
    for j = 1:N
        bases = fBasis(t(j),tf,M);
        x(j) = fourCoeff(1,:)*bases(1,:)';
        y(j) = fourCoeff(2,:)*bases(1,:)';
        z(j) = fourCoeff(3,:)*bases(1,:)';
        
        dx(j) = fourCoeff(1,:)*bases(2,:)';
        dy(j) = fourCoeff(2,:)*bases(2,:)';
        dz(j) = fourCoeff(3,:)*bases(2,:)';
        
        ddx(j) = fourCoeff(1,:)*bases(3,:)';
        ddy(j) = fourCoeff(2,:)*bases(3,:)';
        ddz(j) = fourCoeff(3,:)*bases(3,:)'; 
        
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
    
    % inequality constraints
    c(16*N,1) = -1;
    c((1-1)*N+1:1*N) = z + 0.5*ac.b*abs(sin(mu)); % wing tip clearance
    c((2-1)*N+1:2*N) = -z - 100;
    c((3-1)*N+1:3*N) = x - 500;
    c((4-1)*N+1:4*N) = -x - 500;
    c((5-1)*N+1:5*N) = y - 500;
    c((6-1)*N+1:6*N) = -y - 500;
    c((7-1)*N+1:7*N) = CL - 1.17;
    c((8-1)*N+1:8*N) = -CL - 0.2;
    c((9-1)*N+1:9*N) = CT - 1e-4;
    c((10-1)*N+1:10*N) = -CT - 1e-4;
    c((11-1)*N+1:11*N) = V - 80;
    c((12-1)*N+1:12*N) = -V + 10;
    c((13-1)*N+1:13*N) = mu - pi/3;
    c((14-1)*N+1:14*N) = -mu - pi/3;
    c((15-1)*N+1:15*N) = gam - pi/4;
    c((16-1)*N+1:16*N) = -gam - pi/4;
    
%     % equality constraints
%     ceq(1) = x(1);
%     ceq(2) = y(1);
    ceq = 0;
    
    if nargout > 2 % gradient of the constraints
      dc = [];
      dceq = [];
    end        
    
end