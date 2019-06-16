function [eigE, eigVec, dsSoln, d, N] = getEigStuff(Z, ac)
    N = (length(Z)-3)/16;
    
    VR = Z(end-1);
    tfin = Z(end);
    Psi0 = Z(end-2);
    t = ((tfin/N)*(1:N))';

    % X = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z,df,da,de,dr]
    X(N,16) = 0;
    for j = 1:16
        X(:,j) = Z((j-1)*N+1:j*N);
    end
    X(:,9) = X(:,9) + Psi0*t;

    u = X(:,1); v = X(:,2); w = X(:,3);
    p = X(:,4); q = X(:,5); r = X(:,6);
    Phi = X(:,7); Thet = X(:,8); Psi = X(:,9);
    x = X(:,10); y = X(:,11); z = X(:,12);
    df = X(:,13); da = X(:,14); de = X(:,15); dr = X(:,16);
    
    dsSoln.tfin = tfin;
    dsSoln.VR = VR;
    dsSoln.Psi0 = Psi0;
    dsSoln.t = t;
    dsSoln.u = u; dsSoln.v = v; dsSoln.w = w;
    dsSoln.p = p; dsSoln.q = q; dsSoln.r = r;
    dsSoln.Phi = Phi; dsSoln.Thet = Thet; dsSoln.Psi = Psi;
    dsSoln.x = x; dsSoln.y = y; dsSoln.z = z;
    dsSoln.df = df; dsSoln.da = da; dsSoln.de = de; dsSoln.dr = dr;
    
    d = 9;
    
    trajData.T = dsSoln.tfin; trajData.VR = dsSoln.VR;
    trajData.ac = ac; trajData.N = N;

    trajData.D = fourierdiff(N);

    X = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z];
    U = [df,da,de,dr];

    trajData.X = X; trajData.U = U;

    [~, eigE, ~, ~, eigVec] = spectralMethod(trajData);
end