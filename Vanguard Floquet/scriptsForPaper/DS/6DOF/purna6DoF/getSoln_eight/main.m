clear all
close all
% 02/05/19

load acmod
ac.p_exp = 1; 

% Z = [u(1:N),v(1:N),w(1:N),p(1:N),q(1:N),r(1:N),Phi(1:N),Thet(1:N),...
%      ...Psi(1:N),x(1:N),y(1:N),z(1:N),df(1:N),da(1:N),de(1:N),dr(1:N),CTx(1:N),CTy(1:N),VR,tfin]
% 18*N+2

load('../../test9/results/reshaped/eightLin_reshape_001');
Zinit(end-2) = [];

% load temp_initguess
% Zinit = Z;
% clear Z;

N = (numel(Zinit)-2)/18;

D = fourierDiff(N);

limits = [];
limits(end+1:end+N,:) = repmat([-Inf, Inf],[N,1]);                   % u
limits(end+1:end+N,:) = repmat([-Inf, Inf],[N,1]);                   % v
limits(end+1:end+N,:) = repmat([-Inf, Inf],[N,1]);                   % w
limits(end+1:end+N,:) = repmat([-Inf, Inf],[N,1]);                   % p
limits(end+1:end+N,:) = repmat([-Inf, Inf],[N,1]);                   % q
limits(end+1:end+N,:) = repmat([-Inf, Inf],[N,1]);                   % r
limits(end+1:end+N,:) = repmat([-pi/3, pi/3],[N,1]);                 % Phi
limits(end+1:end+N,:) = repmat([-pi/3, pi/3],[N,1]);                 % Thet
limits(end+1:end+N,:) = repmat([-2*pi, 2*pi],[N,1]);                 % Psi1
limits(end+1:end+N,:) = repmat([-200, 200],[N,1]);                   % x
limits(end+1:end+N,:) = repmat([-200, 200],[N,1]);                   % y
limits(end+1:end+N,:) = repmat([-200, -0.1],[N,1]);                  % z
limits(end+1:end+N,:) = repmat([-pi/4, pi/4],[N,1]);                 % df
limits(end+1:end+N,:) = repmat([-pi/4, pi/4],[N,1]);                 % da
limits(end+1:end+N,:) = repmat([-pi/4, pi/4],[N,1]);                 % de
limits(end+1:end+N,:) = repmat([-pi/4, pi/4],[N,1]);                 % dr
limits(end+1:end+N,:) = repmat([-1, 1],[N,1]);                       % CTx
limits(end+1:end+N,:) = repmat([-1, 1],[N,1]);                       % CTy
limits(end+1,:) = [0, 200];                                          % VR
limits(end+1,:) = [1, 25];                                           % tfin

lb = limits(:,1); ub = limits(:,2);

Aeq(2*N,18*N+2) = 0; beq(2*N,1) = 0;
for i = 1:2*N
    Aeq(i,16*N+i) = 1;
end

% optmization call
% alg = 'interior-point';
alg = 'sqp';
options = optimoptions(@fmincon,'Algorithm',alg,'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,...
                       'Display','iter','UseParallel',true);
options.MaxFunctionEvaluations = 1000000000; % Default: 25000
options.StepTolerance = 1e-6;
options.MaxIterations = 3000; % Default: 1000
[Z,VR,exitflag] = fmincon(@(Z) costfun(Z,N),Zinit,[],[],Aeq,beq,lb,ub,@(Z) Cfun(Z,ac,N,D),options);

%%
save('temp_initguess','Z');

view_results