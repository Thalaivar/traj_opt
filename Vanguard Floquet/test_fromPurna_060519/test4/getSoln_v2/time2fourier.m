clear all

% number of harmonics
M = 50;
load soln2
zCoeff = fitFS(t,z,T,M);

N = 1000;
tt = linspace(0,T,N);

x(N,1) = 0;
y(N,1) = 0;
for i = 1:N
    ff = fBasis(tt(i),T,M);
    x(i) = ff(1,:)*zCoeff(1,:)';
    y(i) = ff(1,:)*zCoeff(2,:)';
end

figure
subplot(1,2,1)
hold on
plot(t,z(:,1),'-b');
plot(tt,x,'--r','LineWidth',1);

subplot(1,2,2)
hold on
plot(t,z(:,2),'-b');
plot(tt,y,'--r','LineWidth',1);

save('solnFourier2','zCoeff','a','b','T','M')