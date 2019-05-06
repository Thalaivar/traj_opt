clear all


% test hypothesis on Chicone (2.28)
d = 2;
N = 4; % must be even
[D, x] = fourierDiff(N);
M = zeros(d*N);
% pi - periodic system
t = x/2;
for i = 1:N
    A = model(t(i));
    M(i,i) = A(1,1); M(i,N+i) = A(1,2);
    M(i+N,i) = A(2,1); M(i+N,i+N) = A(2,2);
end
D = [2*D, zeros(N); zeros(N), 2*D];
[eigVec,eigVals] = eig(D - M);
eigVals = -diag(eigVals);
FTM = get_FTM(t);

FE = log(eig(FTM))/pi;

function A = model(t)
    A = [-1 + 1.5*(cos(t))^2, 1 - 1.5*sin(t)*cos(t);
         -1 - 1.5*sin(t)*cos(t), -1 + 1.5*(sin(t))^2];
end

% to calculate FTM
function K = friedmann_K(h, t)
    A_psi = model(t);
    A_psi_h_by_2 = model(t+0.5*h);
    A_psi_h = model(t+h);
    E = A_psi_h_by_2*(eye(2) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(2) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(2) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(2) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end

function FTM = get_FTM(t)
    h = t(2) - t(1);
    FTM = eye(2);
    for i = 1:length(t)-2
        K = friedmann_K(h, (t(end) - i*h));
        FTM = FTM*K;
    end
end