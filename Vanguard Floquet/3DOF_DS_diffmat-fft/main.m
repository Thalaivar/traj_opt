clear all
close all

load acmod_001
% ac.NCT = 1;
% save('acmod_001','ac');

[Xinit,N,which_diff] = init_guess();

% differentiation matrices
column1 = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*pi/N)]';
D = toeplitz(column1,column1([1, N:-1:2])); % first derivative matrix
column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix

% N: number of points where cineq and ceq are imposed 

% % reusing previous solution; expanding previous solution to larger grid
% load temp_initguess
% N = 50;
% Xinit = expander(X,N);
% clear  X

% bounds
limits(3*N+2,2) = 0;
limits(1:2*N,:) = repmat([-200,200],2*N,1);        % x,y
limits(2*N+1:3*N,:) = repmat([-200,-0.1],N,1);     % z
limits(end-1,:) = [0,150];                         % tf
limits(end,:) = [0,200];                           % VR
lb = limits(:,1); ub = limits(:,2);


% optmization call
alg = 'interior-point';
% alg = 'sqp';
options = optimoptions(@fmincon,'Algorithm',alg,'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,'Display','iter');
options.MaxFunctionEvaluations = 1000000000; % Default: 25000
options.StepTolerance = 1e-6;
options.MaxIterations = 500;% Default: 1000
if which_diff
    [X,VR,exitflag] = fmincon(@costfun,Xinit,[],[],[],[],lb,ub,@(X) Cfun_diffmat(X,ac,N,D,DD),options);
else
    [X,VR,exitflag] = fmincon(@costfun,Xinit,[],[],[],[],lb,ub,@(X) Cfun_fft(X,ac,N),options);
end    
%%
save('temp_initguess','X');

view_results