function xeq = find_equilibrium(xguess)

options = optimoptions('fminunc', 'Algorithm', 'trust-region', 'Display', 'None', 'SpecifyObjectiveGradient',true);
[xeq, fval] = fminunc(@(x) eq_objfun(x),xguess,options);
end





