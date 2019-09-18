global fonttype
global fontsize

load('rawMaterial/result_N400_OLin.mat');
prm.model = '3dof';
prm.p = 1;
prm.correctChi = false; prm.loiter = true;
trajData = getTrajData(Z, 400, prm);
d = 3;

tFE = timeMarchMethodGen(trajData, @jac3DoFTimeInterp, d);
tFE = sort(tFE, 'ComparisonMethod', 'real'); 
fprintf("VSTM DONE!\n")
%%
% NFSRK4 = linspace(100, 10000, 100);
NFSRK4 = [100, 250, linspace(500, 30000, 60)];
FSRK4FEData = zeros(2, length(NFSRK4));
for i = 1:length(NFSRK4)
    [~,FSGrid] = fourierdiff(NFSRK4(i));
    FE = FSRK4(trajData, FSGrid, @jac3DoFTimeInterp, d);
    FE = sort(FE, 'ComparisonMethod', 'real');
%     FSRK4FEData(1, i) = FE(end);
%     FSRK4FEData(2, i) = FE(end-2);
%     FSRK4FEData(3, i) = FE(end-3);
    FSRK4FEData(:,i) = [FE(end); FE(end-2)];
end
fprintf("FSRK4 DONE!\n")
save('rawMaterial\FSRK4FEDataLin.mat', 'FSRK4FEData', 'NFSRK4');
% %%
NPS = linspace(40, 400, 37);
PSFEData = zeros(2, length(NPS));
for i = 1:length(NPS)
    FE = spectralMethodInterp(trajData, @jac3DoF, NPS(i), d);
%     PSFEData(1,i) = FE(1);
%     PSFEData(2,i) = FE(3);
%     PSFEData(3,i) = FE(4);
    PSFEData(:,i) = [FE(1); FE(3)];
end
save('rawMaterial\PSFEDataLin.mat', 'PSFEData', 'NPS');
fprintf("PSDONE!\n")

%%
% load rawMaterial\FSRK4FEDataLin.mat
% load rawMaterial\PSFEDataLin.mat
% 
% FSKR4discrepData = zeros(3, length(NFSRK4));
% for i = 1:length(NFSRK4)
%     FSKR4discrepData(1,i) = abs(tFE(end) - FSRK4FEData(1,i))/abs(tFE(end));
%     FSKR4discrepData(2,i) = abs(tFE(end-2) - FSRK4FEData(2,i))/abs(tFE(end-2));
% %     FSKR4discrepData(3,i) = abs(tFE(end-3) - FSRK4FEData(3,i))/abs(tFE(end-3));
% end
% 
% PSFEdiscrepData = zeros(3, length(NPS));
% for i = 1:length(NPS)
%     PSFEdiscrepData(1,i) = abs(tFE(end) - PSFEData(1,i))/abs(tFE(end));
%     PSFEdiscrepData(2,i) = abs(tFE(end-2) - PSFEData(2,i))/abs(tFE(end-2));
% %     PSFEdiscrepData(3,i) = abs(tFE(end-3) - PSFEData(3,i))/abs(tFE(end-3));
% end
% 
% %%
% figure
% semilogy(NFSRK4, 100*FSKR4discrepData(1,:), 'or');
% hold on
% semilogy(NFSRK4, 100*FSKR4discrepData(2,:), 'sb');
% % semilogy(NFSRK4, 100*FSKR4discrepData(3,:), 'dr');
% xlabel('$N$', 'Interpreter', 'latex');
% ylabel('$\%$ Relative Error', 'FontSize', fontsize, 'FontName', fonttype, 'Interpreter', 'latex');
% % legend({'$\lambda_1 = 0.0338 + 0.2357\mathrm{i}$', '$\lambda_2 = -0.0013$', '$\lambda_3 = -0.0403$'}, 'Interpreter', 'latex');
% set(gca, 'TicklabelInterpreter', 'latex');
% saveas(gcf, 'plots/accuracyFSRK4.eps', 'epsc');
% 
% figure
% semilogy(NPS, 100*PSFEdiscrepData(1,:), 'or');
% hold on
% semilogy(NPS, 100*PSFEdiscrepData(2,:), 'sb');
% % semilogy(NPS, 100*PSFEdiscrepData(3,:), 'dr');
% xlabel('$N$', 'Interpreter', 'latex');
% % ylabel('$\%$ Relative Error', 'FontSize', fontsize, 'FontName', fonttype, 'Interpreter', 'latex');
% legend({'$\lambda_o^1 = -0.0260 + 0.0891\mathrm{i}$', '$\lambda_o^2 = -0.0499$'}, 'Interpreter', 'latex', 'FontSize', 12);
% ylim([1e-8, 1e0])
% set(gca, 'TicklabelInterpreter', 'latex');
% saveas(gcf, 'plots/accuracyPS.eps', 'epsc')
