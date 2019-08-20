global fonttype
global fontsize

load('rawMaterial\loiterExp_001.mat')
d = 10;

prm.model = '6dof';
prm.p = 0.25;
prm.loiter = true; prm.correctChi = false;
trajData = getTrajData(Z, 50, prm);

tFE = timeMarchMethodGen(trajData, @jac6DoFTimeInterp, d);
tFE = sortrows(tFE, 'descend');
% temp = tFE(3); tFE(3) = tFE(4); tFE(4) = temp;

N = linspace(40, 500, 47);
discrepData = zeros(d, length(N));
for i = 1:length(N)
    FE = spectralMethodInterp(trajData, @jac6DoF, N(i), d);
    FE = sortrows(FE, 'descend');
    for j = 1:d
        discrepData(j,i) = abs(FE(j)-tFE(j))/abs(tFE(j));
    end
end

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
grid minor
hold on
plot(N, discrepData(1,:), '--om')
plot(N, discrepData(2,:), '--sb')
plot(N, discrepData(3,:), '--xr')
plot(N, discrepData(4,:), '--dk')
legend({'$\lambda_1$', '$\lambda_2$', '$\lambda_3$', '$\lambda_4$'}, 'Interpreter', 'latex')
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Relative error', 'FontSize', fontsize, 'FontName', fonttype)
saveas(gcf, 'plots\consistencyVDP.jpg', 'jpg')