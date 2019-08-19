global fonttype;
global fontsize;

load('3_300.mat')

fprintf("\n*******************\nPS method on VDP with mu = %.2f for N:\n", mu);
prm.model = 'vanderpol';
prm.mu = mu;
trajData = getTrajData(X, 300, prm);

N = [50, 150, 300];
[FE1, eig1] = spectralMethodInterp(trajData, @jacVDP, N(1), 2);
[FE2, eig2] = spectralMethodInterp(trajData, @jacVDP, N(2), 2);
[FE3, eig3] = spectralMethodInterp(trajData, @jacVDP, N(3), 2);

fprintf("%d, %d, %d\n*******************\n", N(1), N(2), N(3))

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
scatter(real(eig1), imag(eig1), 'xm')
grid minor
xlabel('Re', fontsize, 14, 'FontName', fonttype)
ylabel('Im', fontsize, 14, 'FontName', fonttype')
saveas(gcf, 'plots\diffNVDP1.jpg', 'jpg')

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
scatter(real(eig2), imag(eig2), 'xm')
grid minor
xlabel('Re', fontsize, 14, 'FontName', fonttype)
ylabel('Im', fontsize, 14, 'FontName', fonttype)
saveas(gcf, 'plots\diffNVDP2.jpg', 'jpg')

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
scatter(real(eig3), imag(eig3), 'xm')
grid minor
xlabel('Re', fontsize, 14, 'FontName', fonttype)
ylabel('Im', fontsize, 14, 'FontName', fonttype)
saveas(gcf, 'plots\diffNVDP3.jpg', 'jpg')