clear all
close all

load temp_initguess
load acmod
ac.p_exp = 0.25;

N = (length(Z)-3)/6;

VR = Z(end-1);
tfin = Z(end);
Psi0 = Z(end-2);
t = ((tfin/N)*(1:N))';
fac = 2*pi/tfin;

x = Z(1:N);  
y = Z(N+1:2*N); 
z = Z(2*N+1:3*N);
Phi = Z(3*N+1:4*N);  
Thet = Z(4*N+1:5*N); 
Psi1 = Z(5*N+1:6*N); 
Psi = Psi1 + Psi0*t;

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

Phihat = fft(Phi); dPhihat =1i*[0:N/2-1 0 -N/2+1:-1]'.*Phihat; % 1 is odd
dPhi = fac*real(ifft(dPhihat));
ddPhihat = ((1i*[0:N/2 -N/2+1:-1]').^2).*Phihat; % 2 is even
ddPhi = fac*fac*real(ifft(ddPhihat));

Thethat = fft(Thet); dThethat =1i*[0:N/2-1 0 -N/2+1:-1]'.*Thethat; % 1 is odd
dThet = fac*real(ifft(dThethat));
ddThethat = ((1i*[0:N/2 -N/2+1:-1]').^2).*Thethat; % 2 is even
ddThet = fac*fac*real(ifft(ddThethat));

Psi1hat = fft(Psi1); dPsi1hat =1i*[0:N/2-1 0 -N/2+1:-1]'.*Psi1hat; % 1 is odd
dPsi = fac*real(ifft(dPsi1hat)) + Psi0;
ddPsi1hat = ((1i*[0:N/2 -N/2+1:-1]').^2).*Psi1hat; % 2 is even
ddPsi = fac*fac*real(ifft(ddPsi1hat));

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

pb_2VT = p*ac.b./(2*VT);
qc_2VT = q*ac.c./(2*VT);
rb_2VT = r*ac.b./(2*VT);


figure
plot3(x,-y,-z,'-b');
hold on
strt = plot3(x(1),-y(1),-z(1),'*g'); % start
fin = plot3(x(end),-y(end),-z(end),'or'); % end
xlabel('x');
ylabel('y');
zlabel('z');
view(61,31)
axis equal
legend([strt,fin],'start','end');
title('Trajectory');
grid on
hold off

figure
subplot(2,3,1)
plot(t,CL,'-k');
xlabel('t');
title('C_L');
xlim([0,t(end)]);

subplot(2,3,2)
plot(t,sqrt(u.^2 + v.^2 + w.^2),'-k');
xlabel('t');
title('V_T');
xlim([0,t(end)]);

subplot(2,3,3)
plot(t,-z,'-k');
xlabel('t');
title('Altitude');
xlim([0,t(end)]);

subplot(2,3,4)
hold on
plot(t,pb_2VT,'-k');
plot(t,repmat(0.1,length(t)),'-r');
plot(t,repmat(-0.1,length(t)),'-r');
hold off
xlabel('t');
title('pb/2V_{T}');
xlim([0,t(end)]);

subplot(2,3,5)
hold on
plot(t,qc_2VT,'-k');
plot(t,repmat(0.03,length(t)),'-r');
plot(t,repmat(-0.03,length(t)),'-r');
hold off
xlabel('t');
title('qc/2V_{T}');
xlim([0,t(end)]);

subplot(2,3,6)
hold on
plot(t,rb_2VT,'-k');
plot(t,repmat(0.25,length(t)),'-r');
plot(t,repmat(-0.25,length(t)),'-r');
hold off
xlabel('t');
title('rb/2V_{T}');
xlim([0,t(end)]);

% controls
figure

subplot(2,3,1)
plot(t,CTy,'-k');
xlabel('t');
title('C_{Ty}');
xlim([0,t(end)]);

subplot(2,3,2)
plot(t,df,'-k');
xlabel('t');
title('\delta_f');
xlim([0,t(end)]);

subplot(2,3,3)
plot(t,da,'-k');
xlabel('t');
title('\delta_a');
xlim([0,t(end)]);

subplot(2,3,4)
plot(t,de,'-k');
xlabel('t');
title('\delta_e');
xlim([0,t(end)]);

subplot(2,3,5)
plot(t,dr,'-k');
xlabel('t');
title('\delta_r');
xlim([0,t(end)]);

subplot(2,3,6)
hold on
plot(t,CTx,'-k');
hold off
xlabel('t');
title('C_{Tx}');
xlim([0,t(end)]);

%%%%%%%%%%%%%%%%%%
figure
subplot(2,3,1)
plot(t,CD,'-ro');
xlabel('t');
title('C_D');
xlim([0,t(end)]);

subplot(2,3,2)
plot(t,aoa*180/pi,'-bo');
xlabel('t');
ylabel('deg');
title('\alpha');
xlim([0,t(end)]);

subplot(2,3,4)
plot(t,Phi*180/pi,'-m');
xlabel('t');
ylabel('deg');
title('\phi');
xlim([0,t(end)]);

subplot(2,3,5)
plot(t,Thet*180/pi,'-m');
xlabel('t');
ylabel('deg');
title('\theta');
xlim([0,t(end)]);

subplot(2,3,6)
plot(t,Psi*180/pi,'-m');
xlabel('t');
ylabel('deg');
title('\psi');
xlim([0,t(end)]);

HR = 20;
VR = VR*(HR^ac.p_exp);
fprintf(' VR = %f \n tfin = %f\n',VR,tfin)

% [t',u,v,w,p,q,r,Phi,Thet,Psi,x,y,z]
% [t',df,da,de,dr,CTx,CTy]