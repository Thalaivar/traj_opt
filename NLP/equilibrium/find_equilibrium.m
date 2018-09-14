function xeq = find_equilibrium(xguess, p)

options = optimoptions('fminunc', 'Algorithm', 'trust-region', 'Display', 'None', 'SpecifyObjectiveGradient',true);
[xeq, fval] = fminunc(@(x) eq_objfun(x, p),xguess,options);
end





