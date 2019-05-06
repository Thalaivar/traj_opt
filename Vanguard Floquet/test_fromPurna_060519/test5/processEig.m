clear all

load eigData

t = linspace(0,solnData.T,N+1); t(1) = [];

%% plotting floquet exponents
figure
hold on
p1 = plot(real(eigE),imag(eigE),'xm');
p2 = plot(real(FE),imag(FE),'or','LineWidth',1,'MarkerSize',10);
limSet = 2*max(abs(real(FE)));
xlim([-limSet,limSet]);
title_str1 = horzcat('VDP Floquent Exponents (mu = ',num2str(solnData.mu),')');
title_str2 = horzcat('N = ',num2str(N));
title({title_str1,title_str2});
grid minor
legend([p1,p2],{'Spectral FE','Time-evolved FE'});

% spurious near the center
% spurRe_center = 1.192; % N=100
% spurRe_center = 1.195; % N=50
% idx = find(abs(real(eigE)-spurRe_center)<1e-3);
% eigE(idx)

% good over the vertical line
goodRe_center = 2.383; % N=100
idx = find(abs(real(eigE)-goodRe_center)<1e-3);
eigE(idx)

PE_norm(1,:,:) = PE(1,:,:)./max(max(sqrt(PE(2,:,:).^2 + PE(1,:,:).^2)));
PE_norm(2,:,:) = PE(2,:,:)./max(max(sqrt(PE(2,:,:).^2 + PE(1,:,:).^2)));


%% plotting periodic eigenvectors
ff1 = figure;
pp1 = uipanel('Parent',ff1,'BorderType','none'); 
pp1.Title = horzcat('Periodic Eigenvectors (N = ',num2str(N),', mu = ',num2str(solnData.mu),')'); 
pp1.TitlePosition = 'centertop'; 
pp1.FontSize = 12;
pp1.FontWeight = 'bold';

MM = 1;

subplot(2,2,1,'parent',pp1)
plot(t,real(PE_norm(1,:,idx(MM))),'-r','LineWidth',1);
title('Re(u_1)/max||u||')
xlim([0,t(end)]);

subplot(2,2,2,'parent',pp1)
plot(t,imag(PE_norm(1,:,idx(MM))),'-r','LineWidth',1);
title('Im(u_1)/max||u||')
xlim([0,t(end)]);

subplot(2,2,3,'parent',pp1)
plot(t,real(PE_norm(2,:,idx(MM))),'-r','LineWidth',1);
title('Re(u_2)/max||u||')
xlim([0,t(end)]);

subplot(2,2,4,'parent',pp1)
plot(t,imag(PE_norm(2,:,idx(MM))),'-r','LineWidth',1);
title('Im(u_2)/max||u||')
xlim([0,t(end)]);
