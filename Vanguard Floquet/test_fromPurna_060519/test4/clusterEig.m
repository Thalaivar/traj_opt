clear all

load eigData

Z = [real(eigM),imag(eigM)];
[idx,C] = kmeans(Z,d,'OnlinePhase','on');

figure
hold on
plot(Z(idx==1,1),Z(idx==1,2),'r.','MarkerSize',10);
plot(Z(idx==2,1),Z(idx==2,2),'b.','MarkerSize',10);

plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3) 
plot(real(FM),imag(FM),'c.','MarkerSize',15,'LineWidth',3)  
axis equal

display(FM);
display(C)