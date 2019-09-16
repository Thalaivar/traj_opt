global fonttype
global fontsize

fprintf("Generating eigenvalue plots for vanderpol and glycolytic oscillator...\n")

load('rawMaterial\2_300.mat')
prm.model = 'vanderpol';
prm.mu = mu;
N = 0.5*(length(X)-1);
trajData = getTrajData(X, N, prm);
[vdpFE, vdp_eigE] = spectralMethodGen(trajData, @jacVDP, 2);
lVDP = Louiville(trajData, @jacVDP);
tVDP = timeMarchMethodGen(trajData, @jacVDPTimeInterp, 2); 

load('rawMaterial\glycoRaw.mat')
prm.model = 'glycolytic';
prm.a = a; prm.b = B;
N = 0.5*(length(X)-1);
trajData = getTrajData(X, N, prm);
[glycoFE, glyco_eigE] = spectralMethodGen(trajData, @jacGlycolytic, 2);
lGlyco = Louiville(trajData, @jacGlycolytic);
tGlyco = timeMarchMethodGen(trajData, @jacGlycolyticTimeInterp, 2);
fprintf("\n*******************\nModel parameters\nN = %d\nvan der pol: mu = %.2f\nglycolyctic: a = %.2f; b = %.2f\n*******************\n", N, mu, a, B)

linewt = 1;
markerSZ = 75;
% figure('units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'PaperPositionMode', 'auto');
figure
scatter(real(vdp_eigE), imag(vdp_eigE), 'xm')
hold on
scatter(real(lVDP), imag(lVDP), markerSZ, 'sg', 'LineWidth', linewt);
scatter(real(tVDP), imag(tVDP), 'ob')
% grid minor
xlabel('Re', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Im', 'FontSize', fontsize, 'FontName', fonttype)
set(gca, 'TicklabelInterpreter', 'latex');
saveas(gcf, 'plots\rawResultsVDP.eps', 'epsc')

% figure('units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'PaperPositionMode', 'auto');
figure
scatter(real(glyco_eigE), imag(glyco_eigE), 'xm')
hold on
scatter(real(lGlyco), imag(lGlyco), markerSZ, 'sg', 'LineWidth', linewt);
scatter(real(tGlyco), imag(tGlyco), 'ob');
% grid minor
xlabel('Re', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Im', 'FontSize', fontsize, 'FontName', fonttype)
lgnd = legend('PS', 'Liouville $\lambda_o^{\neq 0}$', 'ASTM');
temp = [lgnd; lgnd.ItemText];
set(temp, 'FontName', 'Times New Roman')
set(temp, 'Interpreter', 'latex');
ylim([-150, 150])
set(gca, 'TicklabelInterpreter', 'latex');
saveas(gcf, 'plots\rawResultsGlyco.eps', 'epsc')