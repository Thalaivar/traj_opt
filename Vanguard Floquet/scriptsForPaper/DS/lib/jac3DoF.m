function A = jac3DoF(x, u, trajData)
    prm.m = 4.5;
    prm.S = 0.473;
    prm.CD0 = 0.0173;
    prm.CD1 = -0.0337;
    prm.CD2 = 0.0517;
    prm.p_exp = trajData.p;
    % evaluate jacobian
    U = [u, trajData.VR];
    Z = x;
    A = JacEval(Z, U, prm);
    if(trajData.p == 0.25)
        A = [A(1:3,1:3), A(1:3,6); A(6,1:3), A(6,6)];
    elseif(trajData.p == 1)
        A = A(1:3,1:3);
    else
        error("Invalid value for p!")
    end
end