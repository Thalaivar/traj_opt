function X0 = inclinedCircleGuess(N, shape, ac)
    [D,fourierGrid] = fourierdiff(N);
    nState = 12; nControl = 4; nOptim = nState + nControl;
    column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
    DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix
    
    if strcmp(shape, 'loiter')
        VR0 = 0.1; T0 = 10;
        diffFac = 2*pi/T0;
        t = (fourierGrid/diffFac)';
        r = 20;
        x = r*(-1+cos(2*pi*t/T0));
        y = -r*sqrt(2)*sin(2*pi*t/T0);
        z = -r*(1-cos(2*pi*t/T0))-0.11; 
    end
    
    % derivatives of position
    dx = diffFac*D*x; ddx = (diffFac^2)*DD*x;
    dy = diffFac*D*y; ddy = (diffFac^2)*DD*y;
    dz = diffFac*D*z; ddz = (diffFac^2)*DD*z;
    
    % get phi, theta, psi time histories
    phi = zeros(N,1); theta = zeros(N,1); psi = zeros(N,1);
    for i = 1:N
        Wx =  VR0*(-z(i)).^ac.p_exp;
        Wxz = (ac.p_exp*VR0)*((-z(i)).^ac.p_exp)./z(i);
        wind = [Wx,Wxz];
        
        Z = evalZ_3DOF([x(i),y(i),z(i)],[dx(i),dy(i),dz(i);ddx(i),ddy(i),ddz(i)],wind,ac);
        psi(i) = Z(2); theta(i) = Z(3); phi(i) = Z(5);
    end
    psi = unwrap(psi);
    
    % euler angle derivatives
    dphi = diffFac*D*phi; ddphi = (diffFac^2)*DD*phi;
    dtheta = diffFac*D*theta; ddtheta = (diffFac^2)*DD*theta;
    dpsi = diffFac*D*psi; ddpsi = (diffFac^2)*DD*psi;
    
    % get non flat states time histories
    u = zeros(N,1); v = zeros(N,1); w = zeros(N,1); 
    p = zeros(N,1); q = zeros(N,1); r = zeros(N,1); 
    df = zeros(N,1); da = zeros(N,1); de = zeros(N,1); dr = zeros(N,1); 
    CTx = zeros(N,1); CTy = zeros(N,1);
    for i = 1:N
        Wx =  VR0*(-z(i)).^ac.p_exp;
        Wxz = (ac.p_exp*VR0)*((-z(i)).^ac.p_exp)./z(i);
        wind = [Wx,Wxz];
        
        Y  = [phi(i),theta(i),psi(i),x(i),y(i),z(i)];
        dY = [dphi(i),dtheta(i),dpsi(i),dx(i),dy(i),dz(i);ddphi(i),ddtheta(i),ddpsi(i),ddx(i),ddy(i),z(i)];
        Z = evalZ_6DOF(Y,dY,wind,ac);
        
        u(i) = Z(1); v(i) = Z(2); w(i) = Z(3);
        p(i) = Z(4); q(i) = Z(5); r(i) = Z(6);
        df(i) = Z(7); da(i) = Z(8); de(i) = Z(9); dr(i) = Z(10);
        CTx(i) = Z(11); CTy(i) = Z(12);
    end
    
    % construct initial decision vector
    X0 = [u;v;w;p;q;r;phi;theta;psi;x;y;z;df;da;de;dr;T0;VR0];
    if strcmp(shape, 'loiter')
        X0(N*nOptim+3) = -2*pi/T0;
    end
end