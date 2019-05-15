% check variation of spurious modes with different parameters
clearvars
clc
addpath('lib')
d = 4;

%% variation with N for linear wind & O shape
load('trajectories\LOTXX_YY_50.mat');
N = [50, 500]; 
lam1 = spectralFE(N(1), d, ac);
lam2 = spectralFE(N(2), d, ac);
FTM = get_FTM(ac, 'time_march');
FTM = [FTM(1:3,1:3),FTM(1:3,6);FTM(6,1:3),FTM(6,6)];
FM = eig(FTM);
subplot(2,2,1)
label1 = 'N = 50'; label2 = 'N = 500';
ptitle = 'Linear wind and O shape';
plot_result(lam1, lam2, label1, label2, ptitle, FM);
%% variation with N for linear wind & 8 shape
load('trajectories\LETXX_YY_50.mat');
N = [50, 500]; 
lam1 = spectralFE(N(1), d, ac);
lam2 = spectralFE(N(2), d, ac);
FTM = get_FTM(ac, 'time_march');
FTM = [FTM(1:3,1:3),FTM(1:3,6);FTM(6,1:3),FTM(6,6)];
FM = eig(FTM);
subplot(2,2,2)
label1 = 'N = 50'; label2 = 'N = 500';
ptitle = 'Linear wind and 8 shape';
plot_result(lam1, lam2, label1, label2, ptitle, FM);
%% variation with N for expo wind & O shape
load('trajectories\EOTXX_YY_50.mat');
N = [50, 500]; 
lam1 = spectralFE(N(1), d, ac);
lam2 = spectralFE(N(2), d, ac);
FTM = get_FTM(ac, 'time_march');
FTM = [FTM(1:3,1:3),FTM(1:3,6);FTM(6,1:3),FTM(6,6)];
FM = eig(FTM);
subplot(2,2,3)
label1 = 'N = 50'; label2 = 'N = 500';
ptitle = 'Exponential wind and O shape';
plot_result(lam1, lam2, label1, label2, ptitle, FM);
%% variation with N for expo wind & 8 shape
load('trajectories\EETXX_YY_50.mat');
N = [50, 500]; 
lam1 = spectralFE(N(1), d, ac);
lam2 = spectralFE(N(2), d, ac);
FTM = get_FTM(ac, 'time_march');
FTM = [FTM(1:3,1:3),FTM(1:3,6);FTM(6,1:3),FTM(6,6)];
FM = eig(FTM);
subplot(2,2,4)
label1 = 'N = 50'; label2 = 'N = 500';
ptitle = 'Exponential wind and 8 shape';
plot_result(lam1, lam2, label1, label2, ptitle, FM);
%% function to calculate spectral FMs
function lam = spectralFE(N,d,ac)
    [D,x] = fourierdiff(N);
    t = x/(2*pi/ac.tf);
    M = zeros(N*d, N*d);
    DD = zeros(N*d, N*d);
    for i = 1:length(t)
        A = ac.get_jac(t(i), 'FD');
        A = [A(1:3,1:3),A(1:3,6);A(6,1:3),A(6,6)];
        for k = 1:d
            kk = (k-1)*N;
            for j = 1:d
                jj = (j-1)*N;
                M(i+kk,i+jj) = A(k,j);
            end
        end
    end
    for i = 1:d
        j = (i-1)*N+1;
        DD(j:j+N-1,j:j+N-1) =2*pi*D/ac.tf;
    end
    lam = -eig(DD-M);
    lam = exp(lam*ac.tf);
end
%% function to plot results
function plot_result(lam1, lam2, label1, label2, ptitle, FM)
    hold on
    s3 = scatter(real(FM), imag(FM), '*r', 'DisplayName', 'True FM', 'LineWidth', 1.5);
    s1 = scatter(real(lam1), imag(lam1), '*m', 'DisplayName', label1);
    s2 = scatter(real(lam2), imag(lam2), 'ob', 'DisplayName', label2);
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), '-k', 'LineWidth', 1, 'DisplayName', 'false');
    title(ptitle);
    grid minor
    axis equal
    legend([s1,s2,s3]);
end