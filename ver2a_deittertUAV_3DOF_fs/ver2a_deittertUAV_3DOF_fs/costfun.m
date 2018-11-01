function [costval,Dcostval] = costfun(X)
    costval = X(end); % VR
    if nargout>1
        Dcostval = [];
    end
end