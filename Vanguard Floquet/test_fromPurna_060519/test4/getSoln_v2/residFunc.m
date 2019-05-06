function residVal = residFunc(Z0,a,b,opts)

T = Z0(end);

[~,z] = ode15s(@(t,z) dynFunc(t,z,a,b),[0,T],Z0(1:2),opts);

residVal = z(end,1) - Z0(1);
residVal(end+1) = z(end,2) - Z0(2);
    
end