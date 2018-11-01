clear all

load acmod_001;
% ac.NCT = 1;
% save('acmod_001','ac');

[Xinit,N,M] = init_guess();

% M: no. of harmonics
% flat outputs: [x,y,z]

% Xinit: (3*(2*M+1)+2)x1
% [fourCoeff of x y z, tf, VR]

% N: number of points where cineq and ceq are imposed 

% load temp_initguess
% Xinit = X;
% clear  X

% bounds
limits(3*(2*M+1)+2,2) = 0;
limits(1:3*(2*M+1),:) = repmat([-500,500],3*(2*M+1),1);
limits(end-1,:) = [0,150];
limits(end,:) = [0,200];
lb = limits(:,1); ub = limits(:,2);

% optimization call
options = optimoptions(@fmincon,'Algorithm','interior-point','SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,'Display','iter');
options.MaxFunctionEvaluations = 1000000000; % Default: 25000
options.StepTolerance = 1e-10;
options.MaxIterations = 10000;% Default: 1000
[X,VR,exitflag] = fmincon(@costfun,Xinit,[],[],[],[],lb,ub,@(X) Cfun(X,ac,M,N),options);

%%
% save('temp_initguess','X');


