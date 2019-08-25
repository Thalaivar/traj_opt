clear all
close all
% 07/06/19

% testing estPhug()

load('solutions/trajectoryOptimized/L50.mat')
addpath('../../lib/')

% setFigProp

% addpath('../floquet_identifier')

[X,U,T,VR] = stateControlMat3DOF(sol, N, 'not-same');
D = fourierdiff(N);
trajData.T = T; trajData.VR = VR;
trajData.D = D; trajData.X = X;
trajData.U = U; trajData.N = N;
trajData.p = p;

[~,eigE,~,~,eigVec] = spectralMethod(trajData);

[grp,ind,trim_ind] = estEigGroups(eigE,N,T,[]);

% diagnostics
% figure
% hold on
% plot(real(eigE),imag(eigE),'xm','LineWidth',1)
% for i = 1:numel(grp)
%     plot(real(grp{i}),imag(grp{i}),'.b','MarkerSize',5);
% end
% x = linspace(1.2*min(real(eigE)),1.2*max(real(eigE)),100);
% plot(x,0.5*pi*N/T*ones(1,100),'-k');
% plot(x,-0.5*pi*N/T*ones(1,100),'-k');
% grid minor
% display(grp)

% eigE(trim_ind) = [];
% eigVec(:,trim_ind) = [];
% 
% phugRef = struct;
% phugRef.mag = 0.5;
% phugRef.phs = 1.593;
% solnRef = struct;
% solnRef.d = 4;
% solnRef.N = N;
% solnRef.T = T;
% solnRef.Vrms = sqrt(sum(X(:,1).^2)/N);
% solnRef.phugRef = phugRef;
% phugMode = estPhug3DoF(eigE,eigVec,ind,solnRef);
% display(phugMode);

rmpath('../../lib')