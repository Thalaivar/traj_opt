addpath('../lib/')
addpath('../')

% load trajectory
load('../trajectories/eight.mat')
NN = N;
T = sol(9*NN+1); VR = sol(9*NN+2);
[~,t] = fourierdiff(NN);
t = T*t/(2*pi);

% interpolate to 100 point grid
N = 100;
X = zeros(9*N,1);
[~,tt] = fourierdiff(N);
tt = T*tt/(2*pi);
for i = 1:9
    j = (i-1)*NN; jj = (i-1)*N;
    X(jj+1:jj+N) = interp_sinc(t, sol(j+1:j+NN), tt);
end

[xdot, xdot_cap] = get_dervs(X, N, T, p, VR);
plot_dervs(xdot, xdot_cap, N, p, T)

rmpath('../lib/')
rmpath('../')

function [xdot, xdot_cap] = get_dervs(X, N, T, p, VR)
    x = zeros(N, 6); 
    u = zeros(N, 3);
    for i = 1:9
        j = (i-1)*N + 1;
        if i <= 6
            x(:,i) = X(j:j+N-1,1);
        else
            u(:,i-6) = X(j:j+N-1,1);
        end
    end
    [D, ~] = fourierdiff(N);
    xdot_cap = zeros(N*6,1);
    for i = 1:6
        j = (i-1)*N+1;
        xdot_cap(j:j+N-1,1) = 2*pi*D*x(:,i)/T;
    end
    
    % derivative from dynamics
    xdot = zeros(N*6,1);
    for i = 1:N
        state = x(i,:)';
        control = u(i,:)';
        xderv = dynamic_model_DS(state, control, p, VR);
        for j = 1:6
            xdot(i+(j-1)*N,1) = xderv(j,1);
        end
    end
end

