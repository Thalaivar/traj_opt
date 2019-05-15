function f = objectiveFunctionFSD(X, type, traj_type)
    if strcmp(traj_type, 'circle')
        N = (length(X)-3)/8;
    else
        N = (length(X)-2)/8;
    end
    if strcmp(type, 'traj')
        f = X(8*N+2,1);
    end
end