global fonttype
global fontsize

load rawMaterial\FEdat.mat

FE_grps = zeros(6,length(NN));
for i=1:length(NN)
    FE = getLam0(grps{i},tfin);
    nn = length(grps{i});

    for j = nn:-1:1
        if real(FE(j))<-10
            FE_grps(1,i) = FE(j);
        end
    end
    
    for j = nn:-1:1
        if real(FE(j))<-2.5 && real(FE(j))>-8
            FE_grps(2,i) = FE(j);
            break
        end
    end
    
    for j = 1:nn
        if real(FE(j))<-0.5 && real(FE(j))>-2
            FE_grps(3,i) = FE(j);
            break
        end
    end
    
    FE_grps(4,i) = FE(3);
    FE_grps(5,i) = FE(2);
    FE_grps(6,i) = FE(1);
    
end
FE_grps = FE_grps.';
FE_grps = FE_grps.';

% idx = ~isnan(real(FE_LMD(1,:)));
idx = real(FE_grps(1,:))<-10;

% perError1 = abs( ( real(FE_LMD(1,idx)) - real(FE_LMD(1,end)) )/real(FE_LMD(1,end)) )*100;
perError1 = abs( FE_grps(1,idx) - FE_grps(1,end) )/abs(FE_grps(1,end))*100;
perError2 = abs( FE_grps(2,:) - FE_grps(2,end) )/abs(FE_grps(2,end))*100;
perError3 = abs( FE_grps(3,:) - FE_grps(3,end) )/abs(FE_grps(3,end))*100;
perError4 = abs( FE_grps(4,:) - FE_grps(4,end) )/abs(FE_grps(4,end))*100;
perError5 = abs( FE_grps(5,:) - FE_grps(5,end) )/abs(FE_grps(5,end))*100;
perError6 = abs( FE_grps(6,:) - FE_grps(6,end) )/abs(FE_grps(6,end))*100;

labels = {};
for i = 1:6
    labels{i} = strcat('$\lambda_', num2str(i), ' = ', num2str(FE_grps(i,end)), '$');
end

% figure
% plot(NN(idx), perError1, '--om')
% grid minor
% xlabel('$N$', 'Interpreter', 'latex', 'FontName', fonttype, 'FontSize', fontsize);
% ylabel('Relative error', 'FontName', fonttype, 'FontSize', fontsize);
% saveas(gcf, 'plots\PSConsistency1.eps', 'eps')
% 
% figure
% plot(NN, perError2, '--dr')
% grid minor
% xlabel('$N$', 'Interpreter', 'latex', 'FontName', fonttype, 'FontSize', fontsize);
% % ylabel('Relative error', 'FontName', fonttype, 'FontSize', fontsize);
% saveas(gcf, 'plots\PSConsistency2.eps', 'eps')
% 
% figure
% plot(NN, perError3, '--sb')
% grid minor
% xlabel('$N$', 'Interpreter', 'latex', 'FontName', fonttype, 'FontSize', fontsize);
% ylabel('Relative error', 'FontName', fonttype, 'FontSize', fontsize);
% saveas(gcf, 'plots\PSConsistency3.eps', 'eps')
% 
% figure
% plot(NN, perError4, '--ok')
% grid minor
% xlabel('$N$', 'Interpreter', 'latex', 'FontName', fonttype, 'FontSize', fontsize);
% % ylabel('Relative error', 'FontName', fonttype, 'FontSize', fontsize);
% saveas(gcf, 'plots\PSConsistency4.eps', 'eps')
% 
% figure
% plot(NN, perError5, '--pc')
% grid minor
% xlabel('$N$', 'Interpreter', 'latex', 'FontName', fonttype, 'FontSize', fontsize);
% ylabel('Relative error', 'FontName', fonttype, 'FontSize', fontsize);
% saveas(gcf, 'plots\PSConsistency5.eps', 'eps')
% 
% figure
% plot(NN, perError6, '--hg')
% grid minor
% xlabel('$N$', 'Interpreter', 'latex', 'FontName', fonttype, 'FontSize', fontsize);
% % ylabel('Relative error', 'FontName', fonttype, 'FontSize', fontsize);
% saveas(gcf, 'plots\PSConsistency6.eps', 'eps')

fprintf("Eigenvalues for N = 400: \n Colour  ->      FE\n")
fprintf("  m  ->  %f + %fi\n", real(FE_grps(1,end)), imag(FE_grps(1,end)))
fprintf("  r  ->  %f + %fi\n", real(FE_grps(2,end)), imag(FE_grps(2,end)))
fprintf("  b  ->  %f + %fi\n", real(FE_grps(3,end)), imag(FE_grps(3,end)))
fprintf("  k  ->  %f + %fi\n", real(FE_grps(4,end)), imag(FE_grps(4,end)))
fprintf("  c  ->  %f + %fi\n", real(FE_grps(5,end)), imag(FE_grps(5,end)))
fprintf("  g  ->  %f + %fi\n", real(FE_grps(6,end)), imag(FE_grps(6,end)))

figure
subplot(3,2,1)
semilogy(NN(idx),perError1,'om','LineWidth',1)
ylim([1e-7,1e0]);
% xlim([190,NN(end)])
ylabel('$\%$ error', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
% xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype);
title(labels{1}, 'Interpreter', 'latex');
% legend(labels{1}, 'Interpreter', 'latex')
% grid minor

subplot(3,2,2)
semilogy(NN,perError2,'or','LineWidth',1)
ylim([1e-7,1e0]);
% xlim([NN(1)-10,NN(end)])
% ylabel('$\%$ error', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
% xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype);
title(labels{2}, 'Interpreter', 'latex');
% grid minor

subplot(3,2,3)
semilogy(NN,perError3,'ob','LineWidth',1)
ylim([1e-7,1e0]);
% xlim([NN(1)-10,NN(end)])
ylabel('$\%$ error', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
% xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype);
title(labels{3}, 'Interpreter', 'latex');
% grid minor

subplot(3,2,4)
semilogy(NN,perError4,'ok','LineWidth',1)
ylim([1e-7,1e0]);
% xlim([NN(1)-10,NN(end)])
% ylabel('$\%$ error', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
% xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype);
title(labels{4}, 'Interpreter', 'latex');
% grid minor

subplot(3,2,5)
semilogy(NN,perError5,'oc','LineWidth',1)
ylim([1e-7,1e0]);
% xlim([NN(1)-10,NN(end)])
ylabel('$\%$ error', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype);
title(labels{5}, 'Interpreter', 'latex');
% grid minor

subplot(3,2,6)
semilogy(NN,perError6,'og','LineWidth',1)
ylim([1e-7,1e0]);
% xlim([NN(1)-10,NN(end)])
% ylabel('$\%$ error', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype)
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', fonttype);
title(labels{6}, 'Interpreter', 'latex');
% grid minor

saveas(gcf, 'plots/consistencyPS6DoF.eps', 'epsc')