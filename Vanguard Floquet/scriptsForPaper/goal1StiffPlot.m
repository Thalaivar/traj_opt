global fonttype
global fontsize
global markerSZ

lineWD = 2; 

load('rawMaterial\4_300.mat')

N = 0.5*(length(X)-1);
fprintf("\n*******************\nGoal 1 : PS performance on stiff systems for VDP \nwith mu = +- %.2f and N = %d\n*******************\n", mu, N);
prm.model = 'vanderpol';
prm.mu = mu;
trajData = getTrajData(X, N, prm);

lFEp = LouivilleVDP(trajData);
[tFEp4, FTMp4] = timeMarchMethodGen(trajData, @sysModel, 2);
[FEp4, eigEp4] = spectralMethodGen(trajData, @jacVDP, 2);
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
scatter(real(eigEp4), imag(eigEp4), 'xm')
hold on
scatter(real(tFEp4), imag(tFEp4), 'ob')
scatter(real(lFEp), imag(lFEp), markerSZ, 'sg', 'LineWidth', lineWD);
grid minor
xlabel('Re', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Im', 'FontSize', fontsize, 'FontName', fonttype)
legend('PS', 'VSTM', 'Liouville')
saveas(gcf, 'plots\VDPplus4.jpg', 'jpg')

load('rawMaterial\m4_300.mat')
prm.mu = mu;
trajData = getTrajData(X, N, prm);

lFEm = LouivilleVDP(trajData);
[tFEm4, FTMm4] = timeMarchMethodGen(trajData, @sysModel, 2);
[FEm4, eigEm4] = spectralMethodGen(trajData, @jacVDP, 2);
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
scatter(real(eigEm4), imag(eigEm4), 'xm')
hold on
scatter(real(tFEm4), imag(tFEm4), 'ob')
scatter(real(lFEm), imag(lFEm), markerSZ, 'sg', 'LineWidth', lineWD);
grid minor
xlabel('Re', 'FontSize', fontsize, 'FontName', fonttype)
ylabel('Im', 'FontSize', fontsize, 'FontName', fonttype)
legend('PS', 'VSTM', 'Liouville')
saveas(gcf, 'plots\VDPminus4.jpg', 'jpg')

function A = sysModel(t, trajData)
    tt = trajData.T*trajData.fourierGrid/(2*pi);
    x = interp_sinc(tt, trajData.X(:,1), t);
    y = interp_sinc(tt, trajData.X(:,2), t);
    A = jacVDP([x,y], trajData);
end

function lam = LouivilleVDP(trajData)
    intQuad = 0;
    for i = 1:trajData.N
        A = jacVDP(trajData.X(i,:), trajData);
        intQuad = intQuad + trace(A);
    end
    
    lam = intQuad/trajData.N;
end