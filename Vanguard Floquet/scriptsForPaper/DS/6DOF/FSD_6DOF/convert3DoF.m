function X0 = convert3DoF(X, N, shape, ac)
    [D,fourierGrid] = fourierdiff(N);
    nState = 12; nControl = 4; nOptim = nState + nControl;
    column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
    DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix
    
    T = X(8*N+1); VR = X(8*N+2);
    diffFac = 2*pi/T;
    t = fourierGrid/diffFac;
    
    if strcmp(shape, 'loiter')
        psiLinearTerm = X(8*N+3);
    end
    
    % 3-D position
    x = X(3*N+1:4*N); y = X(4*N+1:5*N); z = X(5*N+1:6*N); 
    % derivatives of position
    dx = diffFac*D*x; ddx = (diffFac^2)*DD*x;
    dy = diffFac*D*y; ddy = (diffFac^2)*DD*y;
    dz = diffFac*D*z; ddz = (diffFac^2)*DD*z;
    
    % euler angles
    phi = X(7*N+1:8*N); psi = X(N+1:2*N); theta = X(2*N+1:3*N);
    % euler angle derivatives
    dphi = diffFac*D*phi; ddphi = (diffFac^2)*DD*phi;
    dtheta = diffFac*D*theta; ddtheta = (diffFac^2)*DD*theta;
    dpsi = diffFac*D*psi + psiLinearTerm*ones(N,1); ddpsi = (diffFac^2)*DD*psi;
    
     % get non flat states time histories
    u = zeros(N,1); v = zeros(N,1); w = zeros(N,1); 
    p = zeros(N,1); q = zeros(N,1); r = zeros(N,1); 
    df = zeros(N,1); da = zeros(N,1); de = zeros(N,1); dr = zeros(N,1); 
    CTx = zeros(N,1); CTy = zeros(N,1);
    for i = 1:N
        Wx =  VR*(-z(i)).^ac.p_exp;
        Wxz = (ac.p_exp*VR)*((-z(i)).^ac.p_exp)./z(i);
        wind = [Wx,Wxz];
        
        Y  = [phi(i),theta(i),psi(i),x(i),y(i),z(i)];
        dY = [dphi(i),dtheta(i),dpsi(i),dx(i),dy(i),dz(i);ddphi(i),ddtheta(i),ddpsi(i),ddx(i),ddy(i),ddz(i)];
        Z = evalZ_6DOF(Y,dY,wind,ac);
        
        u(i) = Z(1); v(i) = Z(2); w(i) = Z(3);
        p(i) = Z(4); q(i) = Z(5); r(i) = Z(6);
        df(i) = Z(7); da(i) = Z(8); de(i) = Z(9); dr(i) = Z(10);
        CTx(i) = Z(11); CTy(i) = Z(12);
    end
    
    % construct initial decision vector
    X0 = [u;v;w;p;q;r;phi;theta;psi;x;y;z;df;da;de;dr;T;VR];
    if strcmp(shape, 'loiter')
        X0(N*nOptim+3) = psiLinearTerm;
    end
end