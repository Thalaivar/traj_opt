function [FE, grpSize] = spectralMethodModeID(trajData, floqRef)
    % change this parameter according to model
    d = 9;
    % spectral method runs on same number of grid points as trajectory
    D = trajData.D; N = trajData.N; T = trajData.T;
    
    Dmat = zeros(d*N); Mmat = zeros(d*N);
    for i = 1:N
%         A = sysModel(trajData.p, trajData.X(i,:), [trajData.U(i,:), trajData.VR]);
        A = sysModel(trajData.ac, trajData.X(i,:), trajData.U(i,:), trajData.VR);
        for j = 1:d
           for k = 1:d
                Mmat(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
           end
           Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
        end
    end
    
    [eigVec, eigE] = eig(Dmat - Mmat);
    eigE = -1*diag(eigE);
    
    urms = sqrt(sum(trajData.X(:,1).^2)/N);
    
    [FE, grpSize] = modeID(eigE, eigVec, d, N, T, floqRef, urms);
%     [FE,groupSizes,AM] = identifyFloquet(eigE, N, trajData.T);
end

function A = sysModel(ac, Z, U, VR)
    % Z = [u,v,w,p,q,r,phi,theta,psi,x,y,z]
    % U = [df, da, de, dr]
    u = Z(1); v = Z(2); w = Z(3); p = Z(4); q = Z(5); r = Z(6); 
    phi = Z(7); theta = Z(8); psi = Z(9); x = Z(10); y = Z(11); z = Z(12);
    df = U(1); da = U(2); de = U(3); dr = U(4); CTx = 0; CTy = 0;

    Xval = [u,v,w,p,q,r,phi,theta,psi,x,y,z,df,da,de,dr,CTx,CTy].';

    A = numJacEval(ac,Xval,VR); % finite-diff Jacobian

    % d=9
    A(:,10:12) = [];
    A(10:12,:) = [];
end