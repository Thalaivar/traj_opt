function A = jac6DoFTimeInterp(t, trajData)
    ac = trajData.ac;
    VR = trajData.VR;
    
    tt  = trajData.T*trajData.fourierGrid/(2*pi);
    
    % Z = [u,v,w,p,q,r,phi,theta,psi,x,y,z]
    % U = [df, da, de, dr]
    u     = interp_sinc(tt, trajData.X(:,1), t);
    v     = interp_sinc(tt, trajData.X(:,2), t);
    w     = interp_sinc(tt, trajData.X(:,3), t);
    p     = interp_sinc(tt, trajData.X(:,4), t);
    q     = interp_sinc(tt, trajData.X(:,5), t);
    r     = interp_sinc(tt, trajData.X(:,6), t);
    phi   = interp_sinc(tt, trajData.X(:,7), t);
    theta = interp_sinc(tt, trajData.X(:,8), t);
    psi   = interp_sinc(tt, trajData.X(:,9), t);
    x     = interp_sinc(tt, trajData.X(:,10), t);
    y     = interp_sinc(tt, trajData.X(:,11), t);
    z     = interp_sinc(tt, trajData.X(:,12), t);
    df    = interp_sinc(tt, trajData.U(:,1), t);
    da    = interp_sinc(tt, trajData.U(:,2), t);
    de    = interp_sinc(tt, trajData.U(:,3), t);
    dr    = interp_sinc(tt, trajData.U(:,4), t);
    CTx = 0; CTy = 0;
    
    if isfield(trajData, 'chiLinearTerm')
        psi = psi + trajData.chiLinearTerm*t;
    end

    Xval = [u,v,w,p,q,r,phi,theta,psi,x,y,z,df,da,de,dr,CTx,CTy].';

    A = numJacEval(ac,Xval,VR); % finite-diff Jacobian

    if(trajData.p == 0.25)
        A(:,10:11) = [];
        A(10:11,:) = [];
    elseif(trajData.p == 1)
        A(:,10:12) = [];
        A(10:12,:) = [];
    else
        error("Invalid value for p!")
    end
end