function f = optimizeStabilityFSD(X, trajData, windShear)
    N = trajData.N; fGrid = trajData.fourierGrid;
    
    if strcmp(windShear, 'not-same')
        [x,u,T,VR] = stateControlMat3DOF(X,N,windShear);
        trajData.VR = VR;
    elseif strcmp(windShear, 'same')
        [x,u,T] = stateControlMat3DOF(X,N,windShear);
        VR = trajData.VR;
    end
    
%     trajData.X = x; trajData.U = u;
%     trajData.T = T;
%     
%     chiLinearTerm = nan;
%     if strcmp(windShear,'same')
%         if length(X) == 8*N+2, chiLinearTerm = X(8*N+2); end
%     elseif strcmp(windShear, 'not-same')
%         if length(X) == 8*N+3, chiLinearTerm = X(8*N+3); end
%     end
%     
%     if ~isnan(chiLinearTerm)
%         t = T*fGrid/(2*pi);
%         for i = 1:N
%             x(i,2) = x(i,2) + chiLinearTerm*t(i);
%         end
%     end
%     
%     trajData.X = x; trajData.U = u;
    
%     global FE;
    % global eigE;
    % global AM;
    % global groupSizes;
    % global eigVec;

%     [FE, ~, ~, ~, ~] = spectralMethod(trajData);    
%     f = real(FE(1));
    
%    [~,freidmannGrid] = fourierdiff(500);
%    trajData.type = 'full-state';
%     
%    [~,FE] = freidmannMethod(trajData, freidmannGrid);
%    f = max(real(FE));
   f = VR;
end
