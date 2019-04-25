clearvars
clc
load('sailboat_data.mat');
X = [sol(1:N,1), sol(N+1:2*N,1), sol(2*N+1:3*N,1)];
T = sol(3*N+1,1);
%% fit spline to data
t = linspace(0, T, 1000)';
X_spline = zeros(length(t),3);
splinedata.X = X_spline; splinedata.t = t;
for i = 1:3
    X_spline(:,i) = spline(chebt, X(:,i), t);
end
%% setup spectral floquet problem
N = 50; d = 3;

function dy = sysjacobian(t, splinedata)
    
end