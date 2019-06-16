function [c, ceq] = constrainFunctionFSD(X, trajData, windShear)
    p = trajData.p; D = trajData.D; DD = trajData.DD;
    fGrid = trajData.fourierGrid;
    N = trajData.N;
    
    if strcmp(windShear, 'not-same')
        [x,u,T,VR] = stateControlMat3DOF(X,N,windShear);
        trajData.VR = VR;
    elseif strcmp(windShear, 'same')
        [x,u,T] = stateControlMat3DOF(X,N,windShear);
        VR = trajData.VR;
    end
    trajData.T = T;
    
    chiLinearTerm = nan;
    if strcmp(windShear,'same')
        if length(X) == 8*N+2, chiLinearTerm = X(8*N+2); end
    elseif strcmp(windShear, 'not-same')
        if length(X) == 8*N+3, chiLinearTerm = X(8*N+3); end
    end
    
    diffFac = 2*pi/T;
    xdot_cap = zeros(N*6,1);
    for i = 1:6
        j = (i-1)*N+1;
        if i == 2
            if ~isnan(chiLinearTerm)
                xdot_cap(j:j+N-1,1) = diffFac*D*x(:,i) + chiLinearTerm*ones(N,1);
            else
                xdot_cap(j:j+N-1,1) = diffFac*D*x(:,i);
            end 
        else
            xdot_cap(j:j+N-1,1) = diffFac*D*x(:,i);
        end
    end
    
    if ~isnan(chiLinearTerm)
        t = T*fGrid/(2*pi);
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
    if ~isnan(chiLinearTerm)
        % non-periodicity in chi
        ceq(end+1) = chiLinearTerm*T + 2*pi;
    end
    
%     global FE;
%     global eigE;
%     global AM;
%     global groupSizes;
%     global eigVec;
    [FE,~,AM,groupSizes] = spectralMethod(trajData);
    indx = find(imag(FE) == 0);
    
     c(end+1) = FE(indx) + 0.055;
    c(end+1) = -0.06 - FE(indx);
    
    if(indx == 1)
        c(end+1) = 1e-1 - real(FE(indx) - FE(2))/abs(real(FE(indx)));
        c(end+1) = groupSizes(indx) - (0.5*N + 6);
        ceq(end+1) = AM(indx) - 1;
    elseif(indx == 3)
        c(end+1) = 1e-1 - real(FE(1) - FE(indx))/abs(real(FE(indx)));
        c(end+1) = groupSizes(2) - (0.5*N + 6);
        ceq(end+1) = AM(2) - 1;
    else
        error("No complex pair!")
    end

%     if (imag(FE(1)) ~= 0)
%          c(end+1) = 1 - real(FE(1) - FE(3))/abs(real(FE(1)));
%     else
%         c(end+1) = 1 - real(FE(1) - FE(2))/abs(real(FE(1)));
%     end
%     c(end+1) = groupSizes(1) - (N+6);
%     ceq(end+1) = AM(1) - 1;
    
    
%     phugMode = phugoidStuff(eigE,eigVec,N,T,x(:,1));
%     c(end+1) = real(phugMode.domFE) + 0.04;
end

function phugMode = phugoidStuff(eigE, eigVec, N, T, V)
    [~,ind,trim_ind] = estEigGroups(eigE,N,T,[]);
    eigE(trim_ind) = [];
    eigVec(:,trim_ind) = [];
    phugRef = struct;
    phugRef.mag = 0.5;
    phugRef.phs = 1.593;
    solnRef = struct;
    solnRef.d = 3;
    solnRef.N = N;
    solnRef.T = T;
    solnRef.Vrms = sqrt(sum(V.^2)/N);
    solnRef.phugRef = phugRef;
    phugMode = estPhug3DoF(eigE,eigVec,ind,solnRef);
end