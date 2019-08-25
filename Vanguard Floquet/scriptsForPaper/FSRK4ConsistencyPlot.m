global fonttype
global fontsize

load('rawMaterial\result_N400_OLin.mat')
prm.model = '3dof'; prm.loiter = true; prm.correctChi = false;
prm.p = 1;
trajData = getTrajData(Z, N, prm);
d = 3;

[~, FSGrid] = fourierdiff(400);
finalFE = FSRK4(trajData, FSGrid, @jac3DoFTimeInterp, d);

N = linspace(40, 400, 37);
FSRK4FEData = zeros(d, length(N));
for i = 1:length(N)
    [~,FSGrid] = fourierdiff(N(i));
    FE = FSRK4(trajData, FSGrid, @jac3DoFTimeInterp, d);
    FSRK4FEData(:,i) = FE';
end

discrepData = zeros(d, length(N));
for i = 1:length(N)
    discrepData(1,i) = abs(finalFE(1) - FSRK4FEData(2,i))/abs(finalFE(1));
    discrepData(2,i) = abs(finalFE(2) - FSRK4FEData(1,i))/abs(finalFE(2));
    discrepData(3,i) = abs(finalFE(3) - FSRK4FEData(3,i))/abs(finalFE(3));
end

% figure('units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'PaperPositionMode', 'auto');
figure
% grid minor
semilogy(N, 100*discrepData(1,:), '--om');
hold on
semilogy(N, 100*discrepData(2,:), '--sb');
semilogy(N, 100*discrepData(3,:), '--dr');
% legend({'$\lambda_1 = -0.0260 + 0.0891\mathrm{i}$', '$\lambda_2 = -0.0260 - 0.0891\mathrm{i}$', '$\lambda_3 = -0.0499$'}, 'Interpreter', 'latex')
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('$\%$ Relative Error', 'FontSize', fontsize, 'FontName', fonttype, 'Interpreter', 'latex')
saveas(gcf, 'plots\consistencyFSRK4.eps', 'epsc')