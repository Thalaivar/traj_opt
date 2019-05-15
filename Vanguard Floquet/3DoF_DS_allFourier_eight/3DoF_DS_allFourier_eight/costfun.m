function [costval,Dcostval] = costfun(Z)
    costval = Z(end-1); % VR
    if nargout>1
        Dcostval = [];
    end
end