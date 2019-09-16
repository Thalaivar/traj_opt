global fonttype
global fontsize

markerSZ = 75;
lineWD = 1; 

load('rawMaterial\4_300.mat')

N = 0.5*(length(X)-1);
fprintf("\n*******************\nGoal 1 : PS performance on stiff systems for VDP \nwith mu = +- %.2f and N = %d\n*******************\n", mu, N);
prm.model = 'vanderpol';
prm.mu = mu;
trajData = getTrajData(X, N, prm);

lFEp = LouivilleVDP(trajData);
[tFEp4, FTMp4] = timeMarchMethodGen(trajData, @jacVDPTimeInterp, 2);
[FEp4, eigEp4] = spectralMethodGen(trajData, @jacVDP, 2);

% figure('units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'PaperPositionMode', 'auto');
figure 
scatter(real(eigEp4), imag(eigEp4), 'xm')
hold on
scatter(real(lFEp), imag(lFEp), markerSZ, 'sg', 'LineWidth', lineWD);
scatter(real(tFEp4), imag(tFEp4), 'ob')
% grid minor
xlabel('Re', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Im', 'FontSize', fontsize, 'FontName', fonttype)
% legend('PS', 'VSTM', 'Liouville')
set(gca, 'TicklabelInterpreter', 'latex');
saveas(gcf, 'plots\VDPplus4.eps', 'epsc')

load('rawMaterial\m4_300.mat')
prm.mu = mu;
trajData = getTrajData(X, N, prm);

lFEm = LouivilleVDP(trajData);
[tFEm4, FTMm4] = timeMarchMethodGen(trajData, @jacVDPTimeInterp, 2);
[FEm4, eigEm4] = spectralMethodGen(trajData, @jacVDP, 2);

% figure('units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'PaperPositionMode', 'auto');
figure 
scatter(real(eigEm4), imag(eigEm4), 'xm')
hold on
scatter(real(lFEm), imag(lFEm), markerSZ, 'sg', 'LineWidth', lineWD);
scatter(real(tFEm4), imag(tFEm4), 'ob')
% grid minor
xlabel('Re', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Im', 'FontSize', fontsize, 'FontName', fonttype)
lgnd = legend('PS', 'Liouville $\lambda_o^{\neq 0}$', 'ASTM');
temp = [lgnd; lgnd.ItemText];
% set(temp, 'FontSize', 14)
set(temp, 'FontName', 'Times New Roman')
set(temp, 'Interpreter', 'latex')
set(gca, 'TicklabelInterpreter', 'latex');
saveas(gcf, 'plots\VDPminus4.eps', 'epsc')