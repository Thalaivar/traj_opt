% 15/05/19
clear all
close all

load acmod
ac.p_exp = 1;

% Z = [V(1:N), chi(1:N), gam(1:N), x(1:N), y(1:N), z(1:N), CL(1:N), mu(1:N), chi0, VR, tfin]
% 8*N+3

Zinit = init_guess(ac); 
N = (length(Zinit)-3)/8;


% load temp_initguess
% Zinit = Z;
% Zinit = expander(Z,64); N=64;
% clear Z;

[D,~] = fourierDiff(N);

% diagnostics
% Xinit = zeros(N,8);
% for i = 1:8
% Xinit(:,i) = Zinit((i-1)*N+1:i*N);
% end

% bounds
limits(8*N+3,2) = 0;
limits(1:N,:) = repmat([10, 80],[N,1]);             % V  
limits(N+1:2*N,:) = repmat([-2*pi, 2*pi],[N,1]);    % chi
limits(2*N+1:3*N,:) = repmat([-pi/4, pi/4],[N,1]);  % gam
limits(3*N+1:4*N,:) = repmat([-200, 200],[N,1]);    % x
limits(4*N+1:5*N,:) = repmat([-200, 200],[N,1]);    % y
limits(5*N+1:6*N,:) = repmat([-200, -0.001],[N,1]); % z
limits(6*N+1:7*N,:) = repmat([-0.2, 1.17],[N,1]);   % CL
limits(7*N+1:8*N,:) = repmat([-pi/3, pi/3],[N,1]);  % mu
limits(8*N+1,:) = [-Inf, Inf];                      % chi0
limits(8*N+2,:) = [0, 200];                         % VR
limits(8*N+3,:) = [0, 150];                         % tfin

lb = limits(:,1); ub = limits(:,2);

% optmization call
% alg = 'interior-point';
alg = 'sqp';
options = optimoptions(@fmincon,'Algorithm',alg,'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,'Display','iter');
options.MaxFunctionEvaluations = 1000000000; % Default: 25000
options.StepTolerance = 1e-6;
options.MaxIterations = 1000;% Default: 1000
[Z,VR,exitflag] = fmincon(@costfun,Zinit,[],[],[],[],lb,ub,@(Z) Cfun(Z,ac,N,D),options);

%%
save('temp_initguess','Z');

view_results
