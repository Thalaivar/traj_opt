clear all
close all
% 28/05/19

load acmod
ac.p_exp = 0.25;

% Z = [x(1:N), y(1:N), z(1:N), Phi(1:N), Thet(1:N), Psi1(1:N), Psi0, VR, tfin]
% 6*N+3

Zinit = init_guess(ac); 
N = (length(Zinit)-3)/6;

% load temp_initguess
% Zinit = Z;
% Zinit = expander(Z,50); N=50;
% clear Z;

% diagnostics
% Xinit = zeros(N,8);
% for i = 1:6
% Xinit(:,i) = Zinit((i-1)*N+1:i*N);
% end

% bounds
limits = [];
limits(end+1:end+N,:) = repmat([-200, 200],[N,1]);    % x
limits(end+1:end+N,:) = repmat([-200, 200],[N,1]);    % y
limits(end+1:end+N,:) = repmat([-200, -0.1],[N,1]);   % z  
limits(end+1:end+N,:) = repmat([-pi/3, pi/3],[N,1]);  % Phi
limits(end+1:end+N,:) = repmat([-pi/4, pi/4],[N,1]);  % Thet
limits(end+1:end+N,:) = repmat([-2*pi, 2*pi],[N,1]);  % Psi1
limits(end+1,:) = [-Inf, Inf];                        % Psi0
limits(end+1,:) = [0, 200];                           % VR
limits(end+1,:) = [1, 150];                           % tfin
 
lb = limits(:,1); ub = limits(:,2);

% optmization call
% alg = 'interior-point';
alg = 'sqp';
options = optimoptions(@fmincon,'Algorithm',alg,'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,...
                       'Display','iter','UseParallel',false);
options.MaxFunctionEvaluations = 1000000000; % Default: 25000
options.StepTolerance = 1e-6;
options.MaxIterations = 1000; % Default: 1000
[Z,VR,exitflag] = fmincon(@costfun,Zinit,[],[],[],[],lb,ub,@(Z) Cfun(Z,ac,N),options);
%%
save('temp_initguess','Z');

view_results