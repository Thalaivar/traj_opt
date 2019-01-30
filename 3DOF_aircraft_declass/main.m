clearvars 

addpath('solutions');
addpath('floquet');
addpath('trajectory');
addpath('task');
addpath('constraint_funcs');
addpath('aircraft');

load('solutions/trajectory_opt/lin_O.mat')
global eigval;

solution.VR = ac.VR; solution.tf = ac.tf; 
solution.coeffs = ac.coeffs; solution.N = ac.N;
solution.objfun_type = 'floq_new';
solution.constFun_type = 'traj';

% adigator setup
setup.order = 1;
setup.numvar = 3*(2*ac.N+1)+2;
setup.objective = 'objfun';
setup.constraint = 'constFun_traj';
setup.auxdata  = solution;
%adioptions = adigatorOptions('UNROLL', 1, 'OVERWRITE', 1);
%setup.options = adioptions;
adifuncs = adigatorGenFiles4Fmincon(setup);
 
%[ac, sol] = optimize_stability(ac, [sol(1:end-2,1);sol(end,1)], p);
[ac, sol] = optimize_stability(sol, ac.N);

% options = optimoptions('fminunc', 'Display', 'iter-detailed', 'StepTolerance', 1e-15);
% sol = fminunc(@(X) objfun_eig(X, ac), ac.tf, options);

rmpath('aircraft');
rmpath('solutions');
rmpath('floquet');
rmpath('trajectory');
rmpath('task');
rmpath('constraint_funcs');