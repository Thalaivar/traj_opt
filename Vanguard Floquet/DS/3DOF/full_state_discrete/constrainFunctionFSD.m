function [c, ceq] = constrainFunctionFSD(X, trajData, windShear)
    p = trajData.p; D = trajData.D; DD = trajData.DD;
    t = trajData.fourierGrid; type = trajData.shape;
    
    if strcmp(windShear, 'same')
        if strcmp(type, 'circle')
            N = (length(X)-2)/8;
            chiLinearTerm = X(8*N+2);
        else
            N = (length(X)-1)/8;
        end
        T = X(8*N+1);
        VR = trajData.VR;
    else
        if strcmp(type, 'circle')
            N = (length(X)-3)/8;
            chiLinearTerm = X(8*N+3);
        else
            N = (length(X)-2)/8;
        end
        T = X(8*N+1); VR = X(8*N+2);
        trajData.VR = VR;
    end
    
    trajData.T = T;
    x = zeros(N, 6); 
    u = zeros(N, 3);
    for i = 1:8
        j = (i-1)*N + 1;
        if i <= 6
            x(:,i) = X(j:j+N-1,1);
        else
            u(:,i-6) = X(j:j+N-1,1);
        end
    end
    
    diffFac = 2*pi/T;
    xdot_cap = zeros(N*6,1);
    for i = 1:6
        j = (i-1)*N+1;
        if i == 2
            if strcmp(type, 'circle')
                xdot_cap(j:j+N-1,1) = diffFac*D*x(:,i) + chiLinearTerm*ones(N,1);
            else
                xdot_cap(j:j+N-1,1) = diffFac*D*x(:,i);
            end 
        else
            xdot_cap(j:j+N-1,1) = diffFac*D*x(:,i);
        end
    end
    
    
    
    if strcmp(type, 'circle')
        t = T*t/(2*pi);
        for i = 1:N
            x(i,2) = x(i,2) + chiLinearTerm*t(i);
        end
    end
    
    trajData.X = x; trajData.U = u;
    
    xdot = zeros(N*6,1);
    for i = 1:N
        state = x(i,:)';
        control = u(i,:)';
        xderv = dynamic_model_DS(state, control, p, VR);
        for j = 1:6
            xdot(i+(j-1)*N,1) = xderv(j,1);
        end
    end
    
    ceq = xdot - xdot_cap;
    hmin = 0; b = 3;
    % z constraints
    c(1:N,1) = x(:,6) + 0.5*b*abs(sin(u(:,2))) + hmin;
    % mu rate^2 constraint
    c(end+1:end+N,1) = diffFac*diffFac*DD*u(:,2) - 20*pi/180;
    c(end+1:end+N,1) = -20*pi/180 - diffFac*diffFac*DD*u(:,2);
    % CL rate constraint
    c(end+1:end+N,1) = diffFac*D*u(:,1) - 0.4;
    c(end+1:end+N,1) = -0.4 - diffFac*D*u(:,1);
    % CL rate constraint
    c(end+1:end+N,1) = diffFac*diffFac*DD*u(:,1) - 0.4;
    c(end+1:end+N,1) = -0.4 - diffFac*diffFac*DD*u(:,1);
    if strcmp(type, 'circle')
        % non-periodicity in chi
        ceq(end+1) = chiLinearTerm*T + 2*pi;
    end
    
    [FE, ~, AM, groupSizes] = spectralMethod(trajData);
    if (imag(FE(1)) ~= 0)
         c(end+1) = 1e-1 - real(FE(1) - FE(3))/abs(real(FE(1)));
    else
        c(end+1) = 1e-1 - real(FE(1) - FE(2))/abs(real(FE(1)));
    end
    c(end+1) = groupSizes(1) - (N+6);
    ceq(end+1) = AM(1) - 1;
end