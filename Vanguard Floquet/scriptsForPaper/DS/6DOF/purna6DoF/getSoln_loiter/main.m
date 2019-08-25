clear all
close all
% 30/05/19

load acmod
ac.p_exp = 1; 

% Z = [u(1:N),v(1:N),w(1:N),p(1:N),q(1:N),r(1:N),Phi(1:N),Thet(1:N),...
%      ...Psi1(1:N),x(1:N),y(1:N),z(1:N),df(1:N),da(1:N),de(1:N),dr(1:N),CTx(1:N),CTy(1:N),Psi0,VR,tfin]
% 18*N+3

% load('../../purna6DoFDF/results/reshaped/loiterLin_reshape_001.mat')
% load VROptSQPLin.mat
% Z = [dsSoln.u;dsSoln.v;dsSoln.w;dsSoln.p;dsSoln.q;dsSoln.r;dsSoln.Phi;dsSoln.Thet;dsSoln.Psi1;dsSoln.x;dsSoln.y;dsSoln.z;dsSoln.df;dsSoln.da;dsSoln.de;dsSoln.dr;dsSoln.Psi0;dsSoln.VR;dsSoln.T];
% load('results\loiterExp_001.mat')

load VROptSQPExp.mat
% load initGuess.mat
Zinit = Z;
% clear Z;

% N = (numel(Zinit)-3)/18;
% Zinit = [Zinit(1:16*N);Zinit(18*N+1:18*N+3)];
N = (numel(Zinit)-3)/16;

limits = [];
limits(end+1:end+N,:) = repmat([0, 50],[N,1]);                   % u
limits(end+1:end+N,:) = repmat([-50, 50],[N,1]);                   % v
limits(end+1:end+N,:) = repmat([-50, 50],[N,1]);                   % w
limits(end+1:end+N,:) = repmat([-1, 1],[N,1]);                   % p
limits(end+1:end+N,:) = repmat([-1, 1],[N,1]);                   % q
limits(end+1:end+N,:) = repmat([-1, 1],[N,1]);                   % r
limits(end+1:end+N,:) = repmat([-pi/3, pi/3],[N,1]);                 % Phi
limits(end+1:end+N,:) = repmat([-pi/4, pi/4],[N,1]);                 % Thet
limits(end+1:end+N,:) = repmat([-2*pi, 2*pi],[N,1]);                 % Psi1
limits(end+1:end+N,:) = repmat([-200, 200],[N,1]);                   % x
limits(end+1:end+N,:) = repmat([-200, 200],[N,1]);                   % y
limits(end+1:end+N,:) = repmat([-200, -0.1],[N,1]);                  % z
limits(end+1:end+N,:) = repmat([-pi/3, pi/3],[N,1]);                 % df
limits(end+1:end+N,:) = repmat([-pi/3, pi/3],[N,1]);                 % da
limits(end+1:end+N,:) = repmat([-pi/3, pi/3],[N,1]);                 % de
limits(end+1:end+N,:) = repmat([-pi/3, pi/3],[N,1]);                 % dr
% limits(end+1:end+N,:) = repmat([-1, 1],[N,1]);                  % CTx
% limits(end+1:end+N,:) = repmat([-1, 1],[N,1]);                  % CTy
limits(end+1,:) = [-Inf, Inf];                                       % Psi0
VR0 = Z(end-1);
limits(end+1,:) = [1.2*VR0, 200];                                          % VR
% limits(end+1,:) =  [0, 200];
limits(end+1,:) = [0, 35];                                          % tfin
lb = limits(:,1); ub = limits(:,2);

D = fourierDiff(N);
column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix

% linear constraints on control rates
A = zeros(16*N, 16*N+3); b = zeros(16*N,1);
dLim = pi/4; ddLim = pi/60;
for i = 1:4
    j = (i-1)*N; k = j+12*N;
    A(j+1:j+N,k+1:k+N) = (2*pi/10)*D;
    b(j+1:j+N) = dLim;
    A(4*N+j+1:4*N+j+N,k+1:k+N) = -(2*pi/10)*D;
    b(4*N+j+1:4*N+j+N) = dLim;
end
for i = 9:12
    j = (i-1)*N; k = (i-9)*N + 12*N;
    A(j+1:j+N,k+1:k+N) = (2*pi/10)*(2*pi/10)*DD;
    b(j+1:j+N) = ddLim;
    A(4*N+j+1:4*N+j+N,k+1:k+N) = -(2*pi/10)*(2*pi/10)*DD;
    b(4*N+j+1:4*N+j+N) = ddLim;
end

% optmization call
% alg = 'interior-point';
alg = 'sqp';
options = optimoptions(@fmincon,'Algorithm',alg,'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,...
                       'Display','iter','UseParallel',true);
options.MaxFunctionEvaluations = 1000000000; % Default: 25000
options.StepTolerance = 1e-12;
options.MaxIterations = 5000; % Default: 1000
[Z,VR,exitflag] = fmincon(@(Z) costfun(Z,N,D,ac),Zinit,A,b,[],[],lb,ub,@(Z) Cfun(Z,ac,N,D),options);

%%
save('temp_initguess','Z');

% view_results