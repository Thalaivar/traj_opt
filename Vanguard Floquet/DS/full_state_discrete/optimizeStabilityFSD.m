function f = optimizeStabilityFSD(X, trajData, windShear, objType)
    if strcmp(objType, 'stability')
        if strcmp(windShear, 'same')
            if strcmp(trajData.shape, 'circle')
                N = (length(X) - 2)/8;
                trajData.chiLinearTerm = X(8*N+2);
            elseif strcmp(trajData.shape, 'eight')
                N = (length(X) - 1)/8;
            end
        else
            if strcmp(trajData.shape, 'circle')
                N = (length(X) - 3)/8;
                trajData.chiLinearTerm = X(8*N+3);
            elseif strcmp(trajData.shape, 'eight')
                N = (length(X) - 2)/8;
            end
            trajData.VR = X(8*N+2);
        end

        trajData.T = X(8*N+1,1);

        x = zeros(N,6); u = zeros(N,3);
        for i = 1:8
            j = (i-1)*N + 1;
            if i <= 6
                x(:,i) = X(j:j+N-1,1);
            else
                u(:,i-6) = X(j:j+N-1,1);
            end
        end
        
        t = trajData.T*trajData.fourierGrid/(2*pi);
        if strcmp(trajData.shape, 'circle')
            for i = 1:N
                x(i,2) = x(i,2) + trajData.chiLinearTerm*t(i);
            end
        end
        trajData.X = x; trajData.U = u;
        
        [~,f] = spectralMethod(trajData, true);
%         FE = real(freidmannMethod(trajData,x,u));
%         f = max(FE(find(abs(FE) > 1e-8)));
%         f = testObjFun(trajData, x, u, N);
%         f = sum(FE);
%         
%         f = 0;
%         for i = 1:length(FE)
%             f = f + atan(10*FE(i));
%         end
%         
    elseif strcmp(objType, 'VR')
        if strcmp(windShear, 'same')
            error("Wind shear kept constant. Cannot optimize for VR");
        else
            if strcmp(trajData.shape, 'circle')
                N = (length(X) - 3)/8;
            elseif strcmp(trajData.shape, 'eight')
                N = (length(X) - 2)/8;
            end
            f = X(8*N+2);
        end
    end
end
