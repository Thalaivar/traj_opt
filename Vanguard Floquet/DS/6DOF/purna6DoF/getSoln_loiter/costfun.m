function [costval,Dcostval] = costfun(Z,N)
    costval = Z(end-1); % VR
%     costval = Z(end-1) + norm(Z(16*N+1:17*N)) + norm(Z(17*N+1:18*N)); % VR + |CTx| + |CTx|
    if nargout>1
        Dcostval = [];
    end
end