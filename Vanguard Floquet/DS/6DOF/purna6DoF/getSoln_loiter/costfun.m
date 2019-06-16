function [costval,Dcostval] = costfun(Z,N,D,ac)
    costval = Z(end-1); % VR
%     costval = Z(end-1) + norm(Z(16*N+1:17*N)) + norm(Z(17*N+1:18*N)); % VR + |CTx| + |CTx|

%     X = zeros(N,12); U = zeros(N,4);
%     for i = 1:12
%         j = (i-1)*N;
%         X(:,i) = Z(j+1:j+N);
%     end
%     for i = 13:16
%         j = (i-1)*N;
%         U(:,i-12) = Z(j+1:j+N);
%     end
%     trajData.N = N; trajData.X = X;
%     trajData.U = U; trajData.D = D;
%     trajData.ac = ac;
%         
%     trajData.T = Z(16*N+3); trajData.VR = Z(16*N+2);
%     
%     global AM;
%     global groupSizes;
%     global FE;
%     
%     [FE, ~, AM, groupSizes] = spectralMethod(trajData);
%     costval = real(FE(1));

    if nargout>1
        Dcostval = [];
    end
end