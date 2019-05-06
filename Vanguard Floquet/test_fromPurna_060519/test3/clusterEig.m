clear all

load eigData

Z = [real(eigM),imag(eigM)];
[idx,C] = kmeans(Z,d-1,'OnlinePhase','on');

figure
hold on
plot(Z(idx==1,1),Z(idx==1,2),'r.','MarkerSize',10);
plot(Z(idx==2,1),Z(idx==2,2),'b.','MarkerSize',10);
plot(Z(idx==3,1),Z(idx==3,2),'k.','MarkerSize',10);
plot(Z(idx==4,1),Z(idx==4,2),'m.','MarkerSize',10);
plot(Z(idx==5,1),Z(idx==5,2),'g.','MarkerSize',10);

plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
 
% plot(real(FM),imag(FM),'ko',...
%      'MarkerSize',15,'LineWidth',3)  

display(FM);