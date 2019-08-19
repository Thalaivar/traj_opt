function A = jac6DoF(Z, U, trajData)
    ac = trajData.ac;
    VR = trajData.VR;
    % Z = [u,v,w,p,q,r,phi,theta,psi,x,y,z]
    % U = [df, da, de, dr]
    u = Z(1); v = Z(2); w = Z(3); p = Z(4); q = Z(5); r = Z(6); 
    phi = Z(7); theta = Z(8); psi = Z(9); x = Z(10); y = Z(11); z = Z(12);
    df = U(1); da = U(2); de = U(3); dr = U(4); CTx = 0; CTy = 0;

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
