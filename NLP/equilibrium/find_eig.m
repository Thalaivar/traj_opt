function [eigval, eigvec] = find_eig(A)
    N = size(A); N = N(1); p = [N, 1];
    eigval = complex(zeros(N,1),0);
    eigvec = complex(zeros(N,N), 0);
    in_guess = [complex(0,0); complex(ones(N,1),0)];
    options1 = optimoptions('fminunc', 'Display', 'Iter', 'Algorithm', 'quasi-newton', 'StepTolerance', 1e-15);
    options2 = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', 'StepTolerance', 1e-15, 'MaxFunctionEvaluations', 20000);
    [x, fval] = fminunc(@(x) objFun(x, p, A), in_guess, options1);
    eigval(1,1) = x(1,1);
    eigvec(:,1) = x(2:end,1);
    p(2) = eigval(1,1);
    in_guess = [eigval(1,1);eigvec(:,1)];
    [x, fval] = fmincon(@(x) objFun(x, p, A), in_guess, [], [], [], [], [], [], @(x) constFun(x,p), options2);
    eigval(2,1) = x(1,1); eigvec(:,2) = x(2:end,1);
    p(2) = eigval(2,1);
    in_guess = [p(2);eigvec(:,1)];
    [x, fval] = fmincon(@(x) objFun(x, p, A), in_guess, [], [], [], [], [], [], @(x) constFun(x,p), options2);
    eigval(2,1) = x(1,1); eigvec(:,2) = x(2:end,1);
        
    
end