% 14/05/19
clear all
close all


addpath('./JacEval');

addpath('./getSoln_loiter');
load acmod

% load('./getSoln_eight/results/eightExp_struct_001')
% ac.p_exp = 0.25;

load('./getSoln_eight/results/eightLin_struct_001')
ac.p_exp = 1;

% load('./getSoln_loiter/results/loiterExp_struct_001')
% ac.p_exp = 0.25;

% load('./getSoln_loiter/results/loiterLin_struct_002')
% ac.p_exp = 1;

global d
d = 10;
N = 100;

T = dsSoln.T;
[D, tau] = fourierDiff(N);
t = tau*T/(2*pi);

Dmat = zeros(d*N);
Mmat = zeros(d*N);
for i = 1:N
    A = sysModel(t(i),dsSoln,ac);
    for j = 1:d
       for k = 1:d
            Mmat(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
       end
       Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
    end
end

tic
[eigVec,eigE] = eig(Mmat-Dmat);
eigE = diag(eigE);
toc
eigM = exp(eigE*T);

% estimate FTM
tic
FTM = get_FTM(t,dsSoln,ac,0);
FM = eig(FTM);
toc
FE = log(FM)/T;

% eigE(imag(eigE)>0.5*pi*N/T) = [];
% eigE(imag(eigE)<-0.5*pi*N/T) = [];
% [~,idx] = sort(real(eigE),'descend');
% sort_eigE = eigE(idx);
% lagVal = 5;
% for i = 1:d*N-lagVal
%     if ~imModulo(sort_eigE(i:i+lagVal),2*pi/T)
%         maxReFE = real(sort_eigE(i));
%         break
%     end
% end
% display(FE);
% display(maxReFE);

rmpath('./JacEval');
rmpath('./getSoln_loiter');

%% FE

figure
hold on
plot(real(eigE),imag(eigE),'xm','LineWidth',1)
plot(real(FE),imag(FE),'.b','MarkerSize',10);
title('Floquet Exponents');
legend('Spectral FE','Time-march FE');
grid minor

% figure
% hold on
% p1 = plot(real(eigE),imag(eigE),'xm');
% p2 = plot(real(FE),imag(FE),'.b','LineWidth',1,'MarkerSize',10);
% limSet = 1.2*max(abs(real(eigE)));
% xlim([-limSet,limSet]);
% title_str1 = horzcat('6DoF DS Floquent Exponents');
% title_str2 = horzcat('N = ',num2str(N));
% title({title_str1,title_str2});
% grid minor
% legend([p1,p2],{'Spectral FE','Time-evolved FE'});


%%

function A = sysModel(t,dsSoln,ac)

u = sincInterp(dsSoln.t,dsSoln.u,t);
v = sincInterp(dsSoln.t,dsSoln.v,t);
w = sincInterp(dsSoln.t,dsSoln.w,t);
p = sincInterp(dsSoln.t,dsSoln.p,t);
q = sincInterp(dsSoln.t,dsSoln.q,t);
r = sincInterp(dsSoln.t,dsSoln.r,t);
Phi = sincInterp(dsSoln.t,dsSoln.Phi,t);
Thet = sincInterp(dsSoln.t,dsSoln.Thet,t);
Psi = sincInterp(dsSoln.t,dsSoln.Psi1,t) + dsSoln.Psi0*t;
x = sincInterp(dsSoln.t,dsSoln.x,t);
y = sincInterp(dsSoln.t,dsSoln.y,t);
z = sincInterp(dsSoln.t,dsSoln.z,t);
df = sincInterp(dsSoln.t,dsSoln.df,t);
da = sincInterp(dsSoln.t,dsSoln.da,t);
de = sincInterp(dsSoln.t,dsSoln.de,t);
dr = sincInterp(dsSoln.t,dsSoln.dr,t);
CTx = 0;
CTy = 0;


Xval = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z,df,da,de,dr,CTx,CTy].';

A = numJacEval(ac,Xval,dsSoln.VR); % finite-diff Jacobian

% d=10
A(:,10:11) = [];
A(10:11,:) = [];

% d=9
% A(:,10:12) = [];
% A(10:12,:) = [];

end

% to calculate FTM
function FTM = get_FTM(t,dsSoln,ac,flg)
global d
    switch flg
        case 0 % time-evolve method
            Phi=zeros(d);
            Phi0=eye(d);
            options = odeset('RelTol',1e-7,'AbsTol',1e-7);
            linDervs = @(t,X) sysModel(t,dsSoln,ac)*X;
            for i=1:d
                [~,X] = ode15s(@(t,X) linDervs(t,X),[0,t(end)],Phi0(:,i),options);
                Phi(:,i)=X(end,:)';
            end
            FTM = Phi;
        case 1 % Friedmann's method
            h = t(2) - t(1);
            FTM = eye(d);
            for i = 1:length(t)-2
                K = friedmann_K(h, (t(end) - i*h),dsSoln,ac);
                FTM = FTM*K;
            end
    end
 
end

function K = friedmann_K(h,t,dsSoln,ac)
global d
    A_psi = sysModel(t,dsSoln,ac);
    A_psi_h_by_2 = sysModel(t+0.5*h,dsSoln,ac);
    A_psi_h = sysModel(t+h,dsSoln,ac);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end