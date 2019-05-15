clear all
load acmod

% load temp_initguess
load eight_exp_001
ac.p_exp = 0.25;

N = (length(Z)-2)/8;

tfin = Z(end);
t = linspace(tfin/N,tfin,N);

V = Z(1:N);
chi = Z(N+1:2*N);
gam = Z(2*N+1:3*N);
x = Z(3*N+1:4*N);
y = Z(4*N+1:5*N);
z = Z(5*N+1:6*N);
CL = Z(6*N+1:7*N);
mu = Z(7*N+1:8*N);

figure
plot3(x,-y,-z,'-b');
hold on
strt = plot3(x(1),-y(1),-z(1),'*g'); % start
fin = plot3(x(end),-y(end),-z(end),'or'); % end
xlabel('x');
ylabel('y');
zlabel('z');
view(61,31)
legend([strt,fin],'start','end');
title('Trajectory');
grid on
hold off

pltmarker = '-k';
figure
subplot(2,3,1)
plot(t,V,pltmarker);
title('V');

subplot(2,3,2)
plot(t,chi,pltmarker);
title('\chi');

subplot(2,3,3)
plot(t,gam,pltmarker);
title('\gamma');

subplot(2,3,4)
plot(t,CL,pltmarker);
title('C_L');

subplot(2,3,5)
plot(t,mu,pltmarker);
title('\mu');

subplot(2,3,6)
plot(t,zeros(size(t)),pltmarker);
title('C_{T}');