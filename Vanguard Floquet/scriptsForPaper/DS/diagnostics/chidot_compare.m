addpath('../lib')
addpath('..')
load('../LET_200.mat')
x = zeros(N, 6); 
u = zeros(N, 3);
T = sol(9*N+1,1);
VR = sol(9*N+2,1);
for i = 1:9
    j = (i-1)*N + 1;
    if i <= 6
        x(:,i) = sol(j:j+N-1,1);
    else
        u(:,i-6) = sol(j:j+N-1,1);
    end
end
chi = zeros(1,N);
[~,t] = fourierdiff(N);
t = T*t/(2*pi);
for i = 1:N
    xderv = dynamic_model_DS(x(i,:)', u(i,:)', p, VR);
    chi(i) = xderv(2);
end
N = 100;
[~,tt] = fourierdiff(N);
tt = T*tt/(2*pi);
chi_cap = zeros(1,N);
x_cap = zeros(N,6);
u_cap = zeros(N,3);
for i = 1:9
    j = (i-1)*N + 1;
    if i <= 6
        x_cap(:,i) = interp_sinc(t, x(:,i), tt); 
    else
        u_cap(:,i-6) = interp_sinc(t, u(:,i-6), tt);
    end
end
for i = 1:N
    xderv = dynamic_model_DS(x_cap(i,:)', u_cap(i,:)', p, VR);
    chi_cap(i) = xderv(2);
end

rmpath('..')
rmpath('../lib/')