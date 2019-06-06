clear all
close all

load temp_initguess
% load('./results/loiterExp_001');
load acmod
ac.p_exp = 0.25;

N = (length(Z)-3)/18;

VR = Z(end-1);
tfin = Z(end);
Psi0 = Z(end-2);
t = ((tfin/N)*(1:N))';

% X = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z,df,da,de,dr,CTx,CTy]
X(N,18) = 0;
for j = 1:18
    X(:,j) = Z((j-1)*N+1:j*N);
end
X(:,9) = X(:,9) + Psi0*t;

u = X(:,1);
v = X(:,2);
w = X(:,3);
p = X(:,4);
q = X(:,5);
r = X(:,6);
Phi = X(:,7);
Thet = X(:,8);
Psi = X(:,9);
x = X(:,10);
y = X(:,11);
z = X(:,12);
df = X(:,13);
da = X(:,14);
de = X(:,15);
dr = X(:,16);
CTx = X(:,17);
CTy = X(:,18);

VT = sqrt(u.^2 + v.^2 + w.^2);
aoa = atan(w./u);
bet = asin(v./VT);

CL = ( ac.CL0 + ac.CLalf*aoa + ac.CLq*(0.5*q*ac.c./VT) + ac.CLdf*df + ac.CLde*de );
CD = ( ac.CD0 + ac.CD1*CL + ac.CD2*CL.^2 );

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