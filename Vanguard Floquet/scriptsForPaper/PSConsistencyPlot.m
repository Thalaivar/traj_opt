global fonttype
global fontsize

load('rawMaterial\result_N400_OLin.mat')
prm.model = '3dof'; prm.loiter = true; prm.correctChi = true;
prm.p = 1;
trajData = getTrajData(Z, N, prm);
d = 3;

finalFE = spectralMethodGen(trajData, @jac3DoF, d);

prm.correctChi = false;
trajData = getTrajData(Z, N, prm);

N = linspace(40, 400, 37);
PSFEData = zeros(d, length(N));
for i = 1:length(N)
    FE = spectralMethodInterp(trajData, @jac3DoF, N(i), d);
    PSFEData(:,i) = FE';
end

discrepData = zeros(d, length(N));
for i = 1:length(N)
    discrepData(1,i) = abs(finalFE(1) - PSFEData(2,i))/abs(finalFE(1));
    discrepData(2,i) = abs(finalFE(2) - PSFEData(1,i))/abs(finalFE(2));
    discrepData(3,i) = abs(finalFE(3) - PSFEData(3,i))/abs(finalFE(3));
end

figure
% grid minor
semilogy(N, 100*discrepData(1,:), 'ob');
hold on
semilogy(N, 100*discrepData(3,:), 'dr');
legend({'$\lambda_o^1 = -0.0260 + 0.0891\mathrm{i}$', '$\lambda_o^2 = -0.0499$'}, 'Interpreter', 'latex', 'FontSize', 12)
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
% ylabel('$\%$ Relative Error', 'FontSize', fontsize, 'FontName', fonttype, 'Interpreter', 'latex')
ylim([1e-8, 1])
set(gca, 'TicklabelInterpreter', 'latex');
saveas(gcf, 'plots\consistencyPS.eps', 'epsc')