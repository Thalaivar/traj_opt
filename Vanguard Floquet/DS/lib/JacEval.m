function Jac = JacEval(Zval,Uval,prm)
    Z = struct('f',Zval,'dZ',ones(6,1));
    U = Uval;
    y = dac3DOFdyn(Z,U,prm);
    Jac = zeros(y.dZ_size);
    idx = sub2ind(y.dZ_size,y.dZ_location(:,1),y.dZ_location(:,2));
    Jac(idx) = y.dZ;    
end