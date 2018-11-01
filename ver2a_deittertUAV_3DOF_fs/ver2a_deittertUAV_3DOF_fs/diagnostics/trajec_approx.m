clear all
load xyz


% guesses
N = 100;
tfg = 10;
tg  = linspace(0,tfg,N);
xg = 35*(-1+cos(2*pi*tg/tfg));      
yg = -35*sqrt(1)*sin(2*pi*tg/tfg);  
zg = -35*(1-cos(2*pi*tg/tfg))-0.11; 




figure
plot3(x,-y,-z,'-b');
hold on
guess = plot3(xg,-yg,-zg,'-k');
strt = plot3(x(1),-y(1),-z(1),'*g'); % start
fin = plot3(x(end),-y(end),-z(end),'or'); % end
xlabel('x');
ylabel('y');
zlabel('z');
view(61,31)
legend([strt,fin,guess],'start','end','guess');
title('Trajectory');
grid on
hold off