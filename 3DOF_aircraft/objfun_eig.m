 function f = objfun_eig(X, ac)
    global eigval; global solution;
    
    ac.tf = X;
    
    % generate decision vector for inner optimisation
    sol = [ac.coeffs(:,1); ac.coeffs(:,2); ac.coeffs(:,3); ac.VR];
    % get optimal periodic trajectory
    [ac, temp_sol] = optimize_traj(ac, sol, ac.p);
    % get monodromy matrix
    C = get_FTM(ac, 'expo');
    
    % 2-norm
    s = svd(C'*C); f = s(1);
    max_eig = max(eig(C));
    if eigval < max_eig, eigval = max_eig; solution = temp_sol; end       

end