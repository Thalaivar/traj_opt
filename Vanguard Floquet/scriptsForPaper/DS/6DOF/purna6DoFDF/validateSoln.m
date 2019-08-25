clear all

load('./results/loiterExp_001');
load acmod
ac.p_exp = 0.25;

N = (length(Z)-3)/6;

VR = Z(end-1);
tfin = Z(end);
Psi0 = Z(end-2);
t = ((tfin/N)*(1:N))';

x = Z(1:N);  
y = Z(N+1:2*N); 
z = Z(2*N+1:3*N);
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

CTx(N,1) = 0; de(N,1) = 0;
df(N,1) = 0; u(N,1) = 0; v(N,1) = 0; w(N,1) = 0;
p(N,1) = 0; q(N,1) = 0; r(N,1) = 0;
VT(N,1) = 0; bet(N,1) = 0; CL(N,1) = 0;
da(N,1) = 0; dr(N,1) = 0; CTy(N,1) = 0; p(N,1) = 0;
aoa(N,1) = 0;  Q(N,1) = 0; CD(N,1) = 0;

for j = 1:N

    Wx = VR*(-z(j)).^ac.p_exp;
    Wxz = (ac.p_exp*VR)*((-z(j)).^ac.p_exp)./z(j);
    wind = [Wx,Wxz];

    X = evalZ_6DOF([Phi(j),Thet(j),Psi(j),x(j),y(j),z(j)],[dPhi(j),dThet(j),dPsi(j),dx(j),dy(j),dz(j);...
    ddPhi(j),ddThet(j),ddPsi(j),ddx(j),ddy(j),ddz(j)],wind,ac);
    u(j) = X(1);
    v(j) = X(2);
    w(j) = X(3);
    p(j) = X(4);
    q(j) = X(5);
    r(j) = X(6);
    df(j) = X(7);
    da(j) = X(8);
    de(j) = X(9);
    dr(j) = X(10);
    CTx(j) = X(11);
    CTy(j) = X(12);

    VT(j) = sqrt(u(j)^2 + v(j)^2 + w(j)^2);
    aoa(j) = atan(w(j)/u(j));
    bet(j) = asin(v(j)/VT(j));
    
    Q(j) = 0.5*1.225*VT(j)*VT(j);

    CL(j) = ( ac.CL0 + ac.CLalf*aoa(j) + ac.CLq*(0.5*q(j)*ac.c/VT(j)) + ac.CLdf*df(j) + ac.CLde*de(j));        
    CD(j) = ( ac.CD0 + ac.CD1*CL(j) + ac.CD2*CL(j)^2 );

end

% X = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z,df,da,de,dr,CTx,CTy]
X(N,18) = 0;
F(N,12) = 0;
for i = 1:N
    X(i,:) = [u(i),v(i),w(i),p(i),q(i),r(i),Phi(i),Thet(i),Psi(i),x(i),y(i),z(i),df(i),da(i),de(i),dr(i),CTx(i),CTy(i)]';
    F(i,:) = ac.stateDervs(X(i,:),VR);
end

D = fourierDiff(N);
dX = (2*pi/tfin)*D*X(:,1:12); dX(:,9) = dPsi;
residVal = dX - F;
fprintf('max |dX-F| = %s\n',norm(residVal));

Zinit = [u;v;w;p;q;r;Phi;Thet;Psi1;x;y;z;df;da;de;dr;CTx;CTy;Psi0;VR;tfin];
save('./results/reshaped/loiterExp_reshape_001','Zinit');