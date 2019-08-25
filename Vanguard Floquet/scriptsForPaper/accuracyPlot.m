global fonttype
global fontsize
global markerSZ

load('rawMaterial\result_N400_OLin_6DoF.mat')
prm.model = '6dof';
prm.p = 1;
prm.correctChi = false; prm.loiter = true;
trajData = getTrajData(Z, 400, prm);
d = 9;

tFE = timeMarchMethodGen(trajData, @jac6DoFTimeInterp, d);
tFE = sort(tFE, 'ComparisonMethod', 'real'); 
fprintf("VSTM DONE!\n")
%%
% NN = linspace(100, 10000, 100);
% FSRK4FEData = zeros(3, length(NN));
% for i = 1:length(NN)
%     [~,FSGrid] = fourierdiff(NN(i));
%     FE = FSRK4(trajData, FSGrid, @jac6DoFTimeInterp, d);
%     FE = sort(FE, 'ComparisonMethod', 'real');
%     FSRK4FEData(1, i) = FE(end);
%     FSRK4FEData(2, i) = FE(end-2);
%     FSRK4FEData(3, i) = FE(end-3);
% end
% fprintf("FSRK4 DONE!\n")
% %%
% N = linspace(40, 400, 37);
% PSFEData = zeros(3, length(N));
% for i = 1:length(N)
%     FE = spectralMethodInterp(trajData, @jac6DoF, N(i), d);
%     PSFEData(1,i) = FE(1);
%     PSFEData(2,i) = FE(3);
%     PSFEData(3,i) = FE(4);
% end
% fprintf("PSDONE!\n")

%%
load rawMaterial\FSRK4FEData.mat
load rawMaterial\PSFEData.mat

FSKR4discrepData = zeros(3, length(NFSRK4));
for i = 1:length(NFSRK4)
    FSKR4discrepData(1,i) = abs(tFE(end) - FSRK4FEData(1,i))/abs(tFE(end));
    FSKR4discrepData(2,i) = abs(tFE(end-2) - FSRK4FEData(2,i))/abs(tFE(end-2));
    FSKR4discrepData(3,i) = abs(tFE(end-3) - FSRK4FEData(3,i))/abs(tFE(end-3));
end

PSFEdiscrepData = zeros(3, length(NPS));
for i = 1:length(NPS)
    PSFEdiscrepData(1,i) = abs(tFE(end) - PSFEData(1,i))/abs(tFE(end));
    PSFEdiscrepData(2,i) = abs(tFE(end-2) - PSFEData(2,i))/abs(tFE(end-2));
    PSFEdiscrepData(3,i) = abs(tFE(end-3) - PSFEData(3,i))/abs(tFE(end-3));
end

%%
figure
semilogy(NFSRK4, 100*FSKR4discrepData(1,:), '--om');
hold on
semilogy(NFSRK4, 100*FSKR4discrepData(2,:), '--sb');
semilogy(NFSRK4, 100*FSKR4discrepData(3,:), '--dr');
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\%$ Relative Error', 'FontSize', fontsize, 'FontName', fonttype, 'Interpreter', 'latex');
% legend({'$\lambda_1 = 0.0338 + 0.2357\mathrm{i}$', '$\lambda_2 = -0.0013$', '$\lambda_3 = -0.0403$'}, 'Interpreter', 'latex');
saveas(gcf, 'plots/accuracyFSRK4.eps', 'epsc');

figure
semilogy(NPS, 100*PSFEdiscrepData(1,:), '--om');
hold on
semilogy(NPS, 100*PSFEdiscrepData(2,:), '--sb');
semilogy(NPS, 100*PSFEdiscrepData(3,:), '--dr');
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\%$ Relative Error', 'FontSize', fontsize, 'FontName', fonttype, 'Interpreter', 'latex');
legend({'$\lambda_1 = 0.0338 + 0.2357\mathrm{i}$', '$\lambda_2 = -0.0013$', '$\lambda_3 = -0.0403$'}, 'Interpreter', 'latex');
saveas(gcf, 'plots/accuracyPS.eps', 'epsc')