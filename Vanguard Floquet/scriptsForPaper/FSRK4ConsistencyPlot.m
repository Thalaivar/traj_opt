global fonttype
global fontsize

load('rawMaterial\result_N400_OLin.mat')
prm.model = '3dof'; prm.loiter = true; prm.correctChi = false;
prm.p = 1;
trajData = getTrajData(Z, N, prm);
d = 3;

tFE = timeMarchMethodGen(trajData, @jac3DoFTimeInterp, d);

N = linspace(40, 400, 37);
FSRK4FEData = zeros(d, length(N));
for i = 1:length(N)
    [~,FSGrid] = fourierdiff(N(i));
    FE = FSRK4(trajData, FSGrid, @jac3DoFTimeInterp, d);
    FSRK4FEData(:,i) = FE';
end

discrepData = zeros(d, length(N));
for i = 1:length(N)
    discrepData(1,i) = abs(tFE(1) - FSRK4FEData(2,i))/abs(tFE(1));
    discrepData(2,i) = abs(tFE(2) - FSRK4FEData(1,i))/abs(tFE(2));
    discrepData(3,i) = abs(tFE(3) - FSRK4FEData(3,i))/abs(tFE(3));
end

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
grid minor
hold on
plot(N, log(discrepData(1,:)), '--om');
plot(N, log(discrepData(2,:)), '--sb');
plot(N, log(discrepData(3,:)), '--dr');
legend({'$\lambda_1 = -0.0260 + 0.0891\mathrm{i}$', '$\lambda_2 = -0.0260 - 0.0891\mathrm{i}$', '$\lambda_3$ = -0.0499'}, 'Interpreter', 'latex')
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('log Relative Error', 'FontSize', fontsize, 'FontName', fonttype)
saveas(gcf, 'plots\consistencyFSRK4.jpg', 'jpg')