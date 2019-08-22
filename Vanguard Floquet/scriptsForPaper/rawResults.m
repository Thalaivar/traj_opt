global fonttype
global fontsize

fprintf("Generating eigenvalue plots for vanderpol and glycolytic oscillator...\n")

load('rawMaterial\2_300.mat')
prm.model = 'vanderpol';
prm.mu = mu;
N = 0.5*(length(X)-1);
trajData = getTrajData(X, N, prm);
[vdpFE, vdp_eigE] = spectralMethodGen(trajData, @jacVDP, 2);

load('rawMaterial\glycoRaw.mat')
prm.model = 'glycolytic';
prm.a = a; prm.b = B;
N = 0.5*(length(X)-1);
trajData = getTrajData(X, N, prm);
[glycoFE, glyco_eigE] = spectralMethodGen(trajData, @jacGlycolytic, 2);

fprintf("\n*******************\nModel parameters\nN = %d\nvan der pol: mu = %.2f\nglycolyctic: a = %.2f; b = %.2f\n*******************\n", N, mu, a, B)

% figure('units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'PaperPositionMode', 'auto');
figure
scatter(real(vdp_eigE), imag(vdp_eigE), 'xm')
grid minor
xlabel('Re', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Im', 'FontSize', fontsize, 'FontName', fonttype)
saveas(gcf, 'plots\rawResultsVDP.eps', 'eps')

% figure('units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'PaperPositionMode', 'auto');
figure
scatter(real(glyco_eigE), imag(glyco_eigE), 'xm')
grid minor
xlabel('Re', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Im', 'FontSize', fontsize, 'FontName', fonttype)
ylim([-150, 150])
saveas(gcf, 'plots\rawResultsGlyco.eps', 'eps')