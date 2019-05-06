% 01/05/2019
clear all

global d a
a = 0.15;
d = 4;
N =  4;
T = pi;

[D, x] = fourierDiff(N);
t = x*T/(2*pi);

Dmat = zeros(d*N);
Mmat = zeros(d*N);
for i = 1:N
    A = sysModel(t(i));
    for j = 1:d
       for k = 1:d
            Mmat(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
       end
       Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
    end
end

tic
[eigVec,eigE] = eig(Dmat-Mmat);
residVal = norm((Dmat-Mmat)*eigVec-eigVec*eigE);
display(residVal)
eigE = -1*diag(eigE);
eigM = exp(eigE*T);
t1 = toc;

% estimate FTM
tic
FTM = get_FTM(t,0);
FM = eig(FTM);
t2 = toc;

FE = log(FM)/T;
% display(FM)

% % if one of the FE is 0
% intglTrace = 0;
% for i=1:length(t)
% intglTrace = intglTrace + trace(sysModel(t(i)));
% end
% FE2 = [0;intglTrace/numel(t)]; 

% periodic eigenvectors
PE(d,N,N*d) = 0;
for i = 1:N*d
    for j = 1:d
        PE(j,1:N,i) = eigVec((j-1)*N+1:j*N,i);
    end
end

save('eigData','FE','eigE','PE','N','d')

%% FM

% figure
% hold on
% thet = linspace(0,2*pi,100);
% plot(cos(thet),sin(thet),'-k','LineWidth',1);
% axis equal
% title_str = horzcat('Linear Wind , ','N = ',num2str(N));
% title({'Floquet Multipliers',title_str});
% 
% tolVal = min(abs(FM))/5; 
% 
% idx = cell(d);
% ix = [];
% for i = 1:d
%     idx{i} = find(abs(FM(i)-eigM)<tolVal);
%     ix = [ix;idx{i}];
% end
% ix = unique(ix);
% ixx = 1:length(eigM);
% ixx(ix) = [];  
% if isempty(ixx)
%     p0 = plot(0,'m','LineWidth',1);
%     p1 = plot(eigM(ix),'ob','LineWidth',1);
% elseif isempty(ix)
%     p0 = plot(eigM(ixx),'sm','LineWidth',1);
%     p1 = plot(0,'b','LineWidth',1);
% else
%     p0 = plot(eigM(ixx),'sm','LineWidth',1);
%     p1 = plot(eigM(ix),'ob','LineWidth',1);
% end
% p2 = plot(real(FM),imag(FM),'*r','LineWidth',1);
% legend([p1,p0,p2],{horzcat('Accurate FM from spectral ','(',num2str(length(ix)),')'),horzcat('Spurious FM from spectral ','(',num2str(length(ixx)),')'),'FM from time-evolve'});

%% FE

figure
hold on
p1 = plot(real(eigE),imag(eigE),'xm');
p2 = plot(real(FE),imag(FE),'or','LineWidth',1,'MarkerSize',10);
limSet = 1.2*max(abs(real(eigE)));
xlim([-limSet,limSet]);
title_str1 = horzcat('Coupled Oscillators Floquent Exponents');
title_str2 = horzcat('N = ',num2str(N));
title({title_str1,title_str2});
grid minor
legend([p1,p2],{'Spectral FE','Time-evolved FE'});

%% Periodic Eigenvectors

% ff1 = figure;
% pp1 = uipanel('Parent',ff1,'BorderType','none'); 
% pp1.Title = horzcat('Periodic Eigenvectors (N = ',num2str(N),')'); 
% pp1.TitlePosition = 'centertop'; 
% pp1.FontSize = 12;
% pp1.FontWeight = 'bold';
% 
% subplot(2,2,1,'parent',pp1)
% hold on
% for i = 1:d*N
%     plot(t,real(PE(1,:,i)));
% end
% title('Re(u_1)')
% xlim([0,t(end)]);
% 
% subplot(2,2,2,'parent',pp1)
% hold on
% for i = 1:d*N
%     plot(t,imag(PE(1,:,i)));
% end
% title('Im(u_1)')
% xlim([0,t(end)]);
% 
% subplot(2,2,3,'parent',pp1)
% hold on
% for i = 1:d*N
%     plot(t,real(PE(2,:,i)));
% end
% title('Re(u_2)')
% xlim([0,t(end)]);
% 
% subplot(2,2,4,'parent',pp1)
% hold on
% for i = 1:d*N
%     plot(t,imag(PE(2,:,i)));
% end
% title('Im(u_2)')
% xlim([0,t(end)]);

%%

function A = sysModel(t)
global a
    A = [             0, 1,               0, 0;
        -(2+a*cos(2*t)), 0,               1, 0;
                      0, 0,               0, 1;
                      1, 0, -(2+a*cos(2*t)), 0];
end

% to calculate FTM
function FTM = get_FTM(t,flg)
global d
    switch flg
        case 0 % time-evolve method
            Phi=zeros(d);
            Phi0=eye(d);
            options = odeset('RelTol',1e-9,'AbsTol',1e-9);
            linDervs = @(t,X) sysModel(t)*X;
            for i=1:d
                [~,X] = ode15s(@(t,X) linDervs(t,X),[0,t(end)],Phi0(:,i),options);
                Phi(:,i)=X(end,:)';
            end
            FTM = Phi;
        case 1 % Friedmann's method
            h = t(2) - t(1);
            FTM = eye(d);
            for i = 1:length(t)-2
                K = friedmann_K(h, (t(end) - i*h));
                FTM = FTM*K;
            end
    end
 
end

function K = friedmann_K(h,t)
global d
    A_psi = sysModel(t);
    A_psi_h_by_2 = sysModel(t+0.5*h);
    A_psi_h = sysModel(t+h);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end